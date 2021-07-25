#!/usr/bin/env python3
import sys
import argparse
import re
import pysam
import pandas as pd
import numpy as np
import multiprocessing as mp
import scipy.stats as stats
from sklearn.metrics import mutual_info_score, roc_curve, auc
from sklearn import linear_model
from sklearn.model_selection import train_test_split
from functools import partial
from collections import defaultdict

__version__='0.1.3'

base_to_number = {"A":1, "C":2, "G":3, "T":4, "N":5}
number_to_base = dict((v, k) for k, v in base_to_number.items())
ntpairs = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N'}


def binary_search(lista, val):
    """
    Find the index of given value on a list.
    For values that do not exit on the list, it gives
    the index of the inmmediate lower value.
    """
    found = False
    first = 0
    last  = len(lista) - 1
    median = last
    while first <= last and not found:
        median = (first + last)//2
        if lista[median] == val:
            found = True
        else:
            if val < lista[median]:
                last = median - 1
            else:
                first = median + 1
    if lista[median] > val:
        return median - 1
    else:
        return median


def ecdf(x):
    '''Calculate the empirical cumulative density function
       ====================
       Return a cumulative density function based
       on the input data.

    '''
    x = np.array(x)
    x.sort()
    n = len(x)
    y = np.concatenate([[0], np.linspace(1/n, 1, n)])

    def childfunc(sample, x, y, sorted=True):
        if not sorted:
            asort = np.argsort(x)
            x = np.take(x, asort, 0)
            y = np.take(y, asort, 0)
        idx = np.searchsorted(x, sample)
        return y[idx]

    return partial(childfunc, x=x, y=y)


def split_cs_string(cs_string):
    return list(
        zip(
            re.sub('[0-9a-z]', ' ', cs_string).split(),
            re.sub('[:*\-+~]', ' ', cs_string).split()
        )
    )


cslenfuncs = {
    ':': int,
    '*': lambda x: 1,
    '+': lambda x: 0,
    '-': len,
    '~': lambda x: int(re.sub('[a-z]', '', x))
}

def cs_to_df(cs_string, pos):
    cs = split_cs_string(cs_string)
    cslist = list()
    for a, b in cs:
        low = pos
        pos += cslenfuncs[a](b)
        high = pos
        cslist.append([low, high, a, b])
    csdf = pd.DataFrame(
        np.row_stack(cslist),
        columns=['low', 'high', 'ope', 'val']
    )
    csdf.loc[:,'low'] = csdf['low'].astype(int)
    csdf.loc[:,'high'] = csdf['high'].astype(int)
    return csdf

def get_mismatch_from_cs(bam_file, chrom,
                         mapq_threshold = 20,
                         min_allele_count=2,
                         drop_non_spliced_read=True,
                         min_dist_from_splice=4):
    """
    Obtain mismatch coordinates from the cs tags in the bam file
    """
    mm_dict = defaultdict(dict)
    sam = pysam.AlignmentFile(bam_file, 'rb')

    for read in sam.fetch(chrom):
        if read.mapq < mapq_threshold:
                continue
        elif read.is_secondary: ### skip secondary reads
                continue
        # 0 based
        pos = read.reference_start
        cs = cs_to_df(read.get_tag('cs'), pos)
        cs_splice = cs.loc[cs['ope']=='~']
        cs_mismatch = cs.loc[cs['ope']=='*']
        if drop_non_spliced_read:
            if cs_splice.shape[0] == 0:
                continue
        if cs_mismatch.shape[0] > 0:
            if cs_splice.shape[0] > 0:
                select_idx = ~cs_mismatch['low'].map(
                    lambda a: np.logical_and(
                        (cs_splice['low'] - min_dist_from_splice) <= a,
                        (cs_splice['low'] + min_dist_from_splice) >= a
                    ).append(
                        np.logical_and(
                            (cs_splice['high'] - min_dist_from_splice) <= a,
                            (cs_splice['high'] + min_dist_from_splice) >= a
                        )
                    ).any()
                )
            else:
                select_idx = cs_mismatch['low'] > -1
            if select_idx.any():
                for ri, row in cs_mismatch.loc[
                        select_idx, ['low', 'val']
                ].iterrows():
                    pos = int(row['low'])
                    REF, ALT = row['val'].upper()
                    try:
                        mm_dict[(chrom, pos, REF)][ALT] += 1
                    except:
                        mm_dict[(chrom, pos, REF)][ALT] = 1
    arr_chrom = list()
    arr_pos = list()
    arr_ref = list()
    for (chrom, pos, REF), ALT_set in mm_dict.items():
        if len(ALT_set) > 2: ### Removes positions that have more than 2 alleles
            seqs = list(ALT_set.keys())
            values = np.array(list(ALT_set.values()))
            ALT_set = {seqs[i]:values[i] for i in values.argsort()[-2:]}
        if any(AC >= min_allele_count for ALT, AC in ALT_set.items()):
            arr_chrom.append(chrom)
            arr_pos.append(pos)
            arr_ref.append(REF)
    arr = pd.DataFrame(
        {'chromosome':arr_chrom,
         'pos':arr_pos,
         'ref':arr_ref}
    )
    return arr


def filter_mismatch_homopoly(genome_fasta, mismatch_df,
                             chrom, homopoly_length=5):
    genome = pysam.FastaFile(genome_fasta)
    keep_list = []
    for ri, row in mismatch_df.iterrows():
        pos = row['pos']
        left = genome.fetch(
            chrom, pos - homopoly_length, pos
        ).upper()
        right = genome.fetch(
            chrom, pos + 1, pos + homopoly_length + 1
        ).upper()
        not_homo = (len(set(left + right)) > 1)
        keep_list.append(not_homo)
    return mismatch_df.loc[keep_list].reset_index(drop=True).copy()


def filter_mismatch_simple_repeat(repeat_file, mismatch_df, chrom):
    def process_repeat(repeat_file, chrom):
        '''
        Find regions with simple repeats
        '''
        f_sr_region = open(repeat_file)
        repeat_start = list()
        repeat_end = defaultdict(int)
        for line in f_sr_region:
            if line.startswith("#"):
                continue
            line = line.split('\t')
            chromosome, start, end = line[0:3]
            start = int(start)
            end = int(end)
            if chrom == chromosome:
                repeat_start.append(start)
                repeat_end[start] = end
        repeat_start_sorted = sorted(repeat_start)
        return repeat_start_sorted, repeat_end

    sr_starts, sr_ends = process_repeat(repeat_file, chrom)
    idx = list()
    for ri, row in mismatch_df.iterrows():
        pos = int(row['pos'])
        pos_index = binary_search(sr_starts, pos)
        repeat_start = sr_starts[pos_index]
        repeat_end = sr_ends[repeat_start]
        if repeat_start <= pos and pos <= repeat_end:
            idx.append(False)
        else:
            idx.append(True)
    return mismatch_df.loc[idx].reset_index(drop=True).copy()


def mark_mismatch_dbsnp(snp_file, mismatch_df, chrom):
    vcf = pysam.VariantFile(snp_file)
    snp_dict = defaultdict(int)
    for a in vcf.fetch(chrom):
        snp_dict[a.start] = 1
    mismatch_df.loc[:,'snp'] = mismatch_df['pos'].map(snp_dict)
    vcf.close()
    return mismatch_df


def get_gtf_gene_strand(gtf, low, high, padding=500):
    if gtf.shape[0] == 0:
        return None
    strand = gtf['strand'][
        (gtf['low'] - padding >= low) &
        (gtf['high'] + padding <= high)
    ].unique()
    if len(strand) >= 1:
        return strand
    else:
        return None


def get_gtf_exon_strand(gtf, intron, padding=10):
    if gtf.shape[0] == 0:
        return None
    if intron.shape[0] == 0:
        return None
    gtf = gtf.reset_index(drop=True)
    intron = intron.reset_index(drop=True)
    intron_low = intron['low'].astype(int).values
    intron_high = intron['high'].astype(int).values
    gtf_low = gtf['low'].values
    gtf_high = gtf['high'].values
    start_i, start_j = np.where(
        (intron_low[:, None] >= gtf_high - padding) &
        (intron_low[:, None] <= gtf_high + padding)
    )
    end_i, end_j = np.where(
        (intron_high[:, None] >= gtf_low - padding) &
        (intron_high[:, None] <= gtf_low + padding)
    )
    splice_list = []
    if len(start_i) > 0:
        read_gtf_splice_s = pd.DataFrame(
            np.column_stack([
                intron[['low']].values[start_i],
                gtf[['strand', 'high']].values[start_j]
            ]),
            columns=[
                'read_pos', 'gtf_strand', 'gtf_pos'
            ]
        )
        splice_list.append(read_gtf_splice_s)
    if len(end_i) > 0:
        read_gtf_splice_e = pd.DataFrame(
            np.column_stack([
                intron[['high']].values[end_i],
                gtf[['strand', 'low']].values[end_j]
            ]),
            columns=[
                'read_pos', 'gtf_strand', 'gtf_pos'
            ]
        )
        splice_list.append(read_gtf_splice_e)
    if len(splice_list) > 0:
        read_gtf_strand = pd.concat(
            splice_list, axis=0
        )
        strand_count = read_gtf_strand.value_counts(
            'gtf_strand'
        ).sort_values(ascending=False)
        return strand_count.index[strand_count.argmax()]
    else:
        return None


def get_gtf_strand(gtf, low, high, intron,
                   gene_padding=500, exon_padding=10):
    gtf_strand = None
    ## gtf gene strand
    gtf_gene_strand = get_gtf_gene_strand(
        gtf=gtf.loc[gtf['feature']=='gene'],
        low=low,
        high=high,
        padding=gene_padding
    )
    if gtf_gene_strand is not None:
        if len(gtf_gene_strand) == 1:
            gtf_strand = gtf_gene_strand[0]
        else:
            gtf_strand = get_gtf_exon_strand(
                gtf=gtf.loc[gtf['feature']=='exon'],
                intron=intron,
                padding=exon_padding
            )
    return gtf_strand


def get_splice_strand(intron):
    splice_pattern = intron['val'].map(
        lambda x: x[0] + x[1] + x[-2] + x[-1]
    )
    splice_pattern_count = splice_pattern.value_counts()
    strand = None
    if splice_pattern_count.shape[0] >= 1:
        pattern = splice_pattern_count.index[
            splice_pattern_count.argmax()
        ]
        if pattern == 'gtag' or pattern == 'gcag' or pattern == 'atac':
            strand = '+'
        elif pattern == 'ctac' or pattern == 'ctgc' or pattern == 'gtat':
            strand = '-'
    return strand


def get_read_strand(seq_strand, gtf, low, high, intron,
                    gene_padding=500, exon_padding=10):
    ## gtf strand
    gtf_strand = get_gtf_strand(
        gtf=gtf,
        low=low,
        high=high,
        intron=intron,
        gene_padding=gene_padding,
        exon_padding=exon_padding
    )
    ## splice strand
    splice_strand = get_splice_strand(intron=intron)
    # read_strand
    read_strand = seq_strand
    if gtf_strand == splice_strand:
        if gtf_strand is not None:
            read_strand = gtf_strand
    elif seq_strand == gtf_strand:
        read_strand = gtf_strand
    elif seq_strand == splice_strand:
        read_strand = splice_strand
    elif gtf_strand is not None:
        read_strand = gtf_strand
    elif splice_strand is not None:
        read_strand = splice_strand
    return read_strand


def correct_read_strand(bam_file, gtf_file, chrom,
                        mapq_threshold=20,
                        gene_padding=500,
                        exon_padding=10,
                        drop_non_spliced_read=True):
    sam = pysam.AlignmentFile(bam_file, 'rb')
    gtf = pysam.TabixFile(gtf_file)
    read_name_list = []
    seq_strand_list = []
    read_strand_list = []
    for read in sam.fetch(reference=chrom):
        if read.mapq < mapq_threshold:
            continue
        # 0 based
        seq_ref_low = read.reference_start
        seq_ref_high = read.reference_end
        cs = cs_to_df(read.get_tag('cs'), seq_ref_low)
        if drop_non_spliced_read:
            if ~(cs['ope'] == '~').any():
                continue
        gtf_list = list()
        for gtf_entry in gtf.fetch(
                chrom,
                seq_ref_low - gene_padding,
                seq_ref_high + gene_padding,
                parser=pysam.asGTF()):
            gtf_list.append(
                [gtf_entry.feature,
                 gtf_entry.gene_name,
                 gtf_entry.start, gtf_entry.end,
                 gtf_entry.strand]
            )
        gtf_df = pd.DataFrame.from_records(
            gtf_list,
            columns=['feature', 'gene_name',
                     'low', 'high', 'strand']
        )
        del gtf_list
        # correct strand
        seq_strand = ['+', '-'][int(read.is_reverse)]
        # read strand
        read_strand = get_read_strand(
            seq_strand, gtf=gtf_df,
            low=seq_ref_low,
            high=seq_ref_high,
            intron=cs.loc[cs['ope'] == '~'],
            gene_padding=gene_padding,
            exon_padding=exon_padding
        )
        read_name_list.append(read.query_name)
        seq_strand_list.append(seq_strand)
        read_strand_list.append(read_strand)
    df = pd.DataFrame(
        {'read_name': read_name_list,
         'seq_strand': seq_strand_list,
         'read_strand': read_strand_list}
    )
    sam.close()
    gtf.close()
    return df


def process_cigar(read):
    """
    Obtained the matched sequence of the read
    create a dictionary for easy conversion of genomic position to read position
    """
    j = 0
    # Process Sequence
    proc_seq  = ''
    # Process Coordinates
    gPos = read.pos
    ConvPos = dict()
    for c, d in read.cigar:
        ConvPos[gPos] = j
        if c == 0: # M matches
            proc_seq  += read.query_sequence[j : j + d]
            gPos += d
            j += d
        elif c in [2, 3]: # D, N deletions
            gPos += d
            continue
        elif c in [1, 4, 5]: # I, S, H clipping
            j += d
        else:
            raise ValueError("Cannot process CIGAR string: {}".format(cigar_dict[c]))
    return proc_seq, ConvPos


def find_read_clusters(bam_file, chrom):
    """
    Find read clusters in the bam file
    """
    Read_starts_and_ends = defaultdict(int)

    bam_handle = pysam.Samfile(bam_file, 'rb')

    for read in bam_handle.fetch(reference = chrom, until_eof = True):
        Read_starts_and_ends[read.pos] += 1
        Read_starts_and_ends[read.positions[-1]] -= 1

    bam_handle.close()

    Cum_Cov = 0
    Max_Cov = 0
    RC_started = False
    RC_coordinates = list()

    for coord in sorted(Read_starts_and_ends):
        Cum_Cov += Read_starts_and_ends[coord]
        if Cum_Cov > Max_Cov:
            Max_Cov = Cum_Cov
        if Cum_Cov >= 1 and not RC_started:
            start_coord = coord
            RC_started = True
        elif Cum_Cov < 1 and RC_started:
            RC_coordinates.append((start_coord, coord, Max_Cov))
            RC_started = False
            Max_Cov = 0

    Read_Clusters = list()
    i = 0
    for (start, end, Max_Cov) in RC_coordinates:
        if Max_Cov >= 2: ### Skipping single-read read clusters
            Read_Clusters.append((i, start, end, Max_Cov))
            i += 1

    del RC_coordinates
    del Read_starts_and_ends
    return Read_Clusters


def get_mm_info_from_bam(bam_file, chrom, start, end, mm_coords, corrected_strand):
    """
    Obtain overlapping reads for the mismatches
    Only select reads that overlap the read cluster
    mismatches that don't overlap the read cluster are skipped
    """
    mm_info = defaultdict(dict)
    ### Parsing mm coordinates
    FirstSearchIndex = 0
    LastSearchIndex = mm_coords.shape[0] - 1
    mm_sorted_position, mm_sorted_ref_base, mm_sorted_issnp = np.transpose(mm_coords)
    mm_sorted_position = mm_sorted_position.astype(int)
    ### Read bam file
    bam_handle  = pysam.Samfile(bam_file, 'rb')

    for read in bam_handle.fetch(reference = chrom, start = start, end = end):
        read_name = read.query_name
        # corrected_strand
        if read_name not in corrected_strand:
            continue
        read_strand = corrected_strand[read_name]

        RL = len(read.query_sequence)

        ### Filter reads

        if read.is_supplementary:
            continue
        if read.is_secondary:
            continue
        if read.mapq < 20:
            continue
        if 'NH' in dict(read.tags) and read.get_tag('NH') > 1:
            continue

        ### Move mismatch index until it overlaps the read

        while FirstSearchIndex <= LastSearchIndex and mm_sorted_position[FirstSearchIndex] < read.pos:
            FirstSearchIndex += 1

        ### Get the first index for mismatch the overlaps the read
        FirstSearchIndexRead = FirstSearchIndex

        ### Process sequence from CIGAR string
        ProcReadSeq, ConvPos = process_cigar(read)

        read_pos = 0

        for sB, eB in read.get_blocks():
            bL = eB - sB
            ### Get the the first mismatch index that overlaps each read block
            while FirstSearchIndexRead <= LastSearchIndex and mm_sorted_position[FirstSearchIndexRead] <= eB - 1:
                if mm_sorted_position[FirstSearchIndexRead] < sB:
                    FirstSearchIndexRead += 1
                    continue

                SNV_POS = int(mm_sorted_position[FirstSearchIndexRead])
                snv_ref_seq = number_to_base[mm_sorted_ref_base[FirstSearchIndexRead]]
                SNV_REF = snv_ref_seq if read_strand == '+' else ntpairs[snv_ref_seq]
                SNV_ISSNP = mm_sorted_issnp[FirstSearchIndexRead]

                mm_read_index = read_pos + SNV_POS - sB
                ### Obtain read position of the SNV
                ### Actually shortest distance to either end of the read
                SNV_READ_POS = ConvPos[sB] + SNV_POS - sB + 1
                if read_strand == '-':
                    SNV_READ_POS = RL - SNV_READ_POS - 1

                ### Obtain read base of the SNV
                seq_strand = '-' if read.is_reverse else '+'
                if seq_strand == read_strand:
                    snv_read_seq = ProcReadSeq[mm_read_index]
                else:
                    snv_read_seq = ntpairs[ProcReadSeq[mm_read_index]]
                SNV_ALT = snv_read_seq if read_strand == '+' else ntpairs[snv_read_seq]

                FirstSearchIndexRead += 1

                ### Dictionary with read-level information
                ############################################
                ### For short-read we collect read position and read direction
                ### May not be relevant feature in long-read
                # using corrected_strand

                mm_info[
                    (SNV_POS, SNV_REF, read_strand, SNV_ISSNP)
                ][read_name] = SNV_ALT

            read_pos += bL

    bam_handle.close()
    return mm_info


def variant_filter(mm_info, chrom,
                   min_AB=0.1, min_AC=3,
                   min_het_snp_ratio=0.35, max_het_snp_ratio=0.65):
    """
    Filtering of alleles
    """
    mm_info_result = defaultdict(dict)
    ratio_dict = defaultdict(float)
    items_list = list(mm_info.items())

    for ((POS, REF, STRAND, ISSNP), read_dict) in items_list:
        allele_counts = dict((base, []) for base in ['A', 'T', 'C', 'G', 'N'])
        for read_name, read_snv in read_dict.items():
                allele_counts[read_snv].append(read_name)

        DP = sum(len(v) for v in allele_counts.values())

        the_ratio = defaultdict(list)
        if ISSNP == 1:
            mismatch_type = 'snp'
        else:
            mismatch_type = 'mismatch'
        for base, read_list in list(allele_counts.items()):
            AC = len(read_list)
            AB = AC/DP
            ### For all SNV positions, discard alleles that have low allele freq.
            ### or low allelic read count by removing the reads that contain these alleles
            if AB < min_AB or AC < min_AC:
                for read in read_list:
                    mm_info[(POS, REF, STRAND, ISSNP)].pop(read)
                allele_counts.pop(base)
            if mismatch_type == 'snp':
                if (AB >= min_het_snp_ratio) and (AB <= max_het_snp_ratio):
                    mismatch_type = 'het_snp'
            the_ratio[base] = [AC, DP, AB]
        ### If only one allele is left at this position (including the REF)
        ### then discard this position.
        ### For the retained SNVs, change dictionary key format
        if len(allele_counts) < 2:
            mm_info.pop((POS, REF, STRAND, ISSNP))
        else:
            coord = mm_info.pop((POS, REF, STRAND, ISSNP))
            ALTS = [base for base in allele_counts if base != REF]
            for ALT in ALTS:
                variant_id = "{}:{}:{}:{}:{}>{}".format(
                    mismatch_type, chrom, POS, STRAND, REF, ALT
                )
                mm_info_result[variant_id] = coord
                ratio_dict[variant_id] = the_ratio[ALT]
    return mm_info_result, ratio_dict


def mutual_information(mm_info, rc_i,
                       mi_testable_common_reads=5,
                       mi_testable_mono=1):

    """
    Calculate the mutual information between a pair of SNVs
    """
    mi_info = dict()

    sorted_SNV_list = sorted(mm_info.keys(), key = lambda x: x[0])

    for snv1 in sorted_SNV_list: ## SNVs 0
        var1_reads = mm_info[snv1]

        MI_values = []
        MI_coverage = []
        SNV_pair_jakarta = []
        SNV_pair = []

        for snv2 in sorted_SNV_list: ## SNPs 1
            if snv2.find('het_snp') < 0:
                ## skip non het_snp
                continue
            var2_reads = mm_info[snv2]

            if snv1 == snv2:
                continue
            common_reads = set(var1_reads) & set(var2_reads)
            lc = len(common_reads)

            ### check overlap of reads between SNVs (Jakarta index)
            jakarta = lc/len(set(var1_reads) | set(var2_reads))


            if lc < mi_testable_common_reads:
                continue

            base_to_number = {"A":1, "C":2, "G":3, "T":4, "N":5}
            mat = np.zeros((lc, 2))
            for i, r in enumerate(common_reads):
                mat[i][0] = base_to_number[var1_reads[r]]
                mat[i][1] = base_to_number[var2_reads[r]]

            var1_alleles, var1_allele_counts = np.unique(mat[:, 0], return_counts = True)
            var2_alleles, var2_allele_counts = np.unique(mat[:, 1], return_counts = True)

            ### Only selecting pairs of SNVs that have mono >= mi_testable_mono (1)
            if any(AC > lc - mi_testable_mono for AC in var1_allele_counts):
                continue

            if any(AC > lc - mi_testable_mono for AC in var2_allele_counts):
                continue

            SNV_pair.append(snv2)
            SNV_pair_jakarta.append(jakarta)

            MI = mutual_info_score(mat[:, 0], mat[:, 1])
            MI_values.append(round(MI, 3))
            MI_coverage.append(lc)

        ### save information of SNVs that were testable for MI
        if MI_values:
            ### weighed MI average by common read count
            wAve_MI  = np.average(MI_values, weights = MI_coverage)
            n = len(MI_coverage)
            mi_info[snv1] = [n, wAve_MI, SNV_pair, SNV_pair_jakarta, MI_values, MI_coverage]

    return mi_info


def format_mi_output_file(mi_info):
    """
    Format output for nice printing
    """
    output_lines = []
    for snv, (n, mi, snv_pairs, snv_pair_jk, mi_vars, mi_cov) in mi_info.items():
        mi_vars     = ["{:.3f}".format(_mi) for _mi in mi_vars]
        snv_pair_jk = ["{:.3f}".format(_jk) for _jk in snv_pair_jk]
        mi = round(mi, 3)
        line = snv.split(':') + [mi, n, ";".join(snv_pairs), ";".join(snv_pair_jk), ";".join(mi_vars), ";".join(map(str, mi_cov))]
        output_lines.append(line)
    return output_lines


def diff_of_allelic_ratio(mm_ratio_df, gtf_file, chrom, default_ratio = 0.5):

    """
    Calculate the difference between editing ratio and gene allelic ratio (mean heterozygous SNP ratio)
    """
    gtf = pysam.TabixFile(gtf_file)

    mm_ratio_df = mm_ratio_df.loc[
        mm_ratio_df['chromosome'] == chrom
    ]
    mm_ratio_df.loc[:, 'mmidx'] = mm_ratio_df.apply(
        lambda x: '{}:{}:{}:{}:{}'.format(
            x['type'], x['chromosome'], x['pos'],
            x['strand'], x['change_type']
        ),
        axis = 1
    )
    mm_gene = defaultdict(lambda : "Not_found")
    gene_het_snp_ratio = defaultdict(list)

    for ri, row in mm_ratio_df.iterrows():
        gtf_list = list()
        for gtf_entry in gtf.fetch(
                chrom,
                row['pos'] - 10,
                row['pos'] + 10,
                parser=pysam.asGTF()):
            gtf_list.append(
                [gtf_entry.feature,
                 gtf_entry.gene_name,
                 gtf_entry.start, gtf_entry.end,
                 gtf_entry.strand]
            )
        if len(gtf_list) > 0:
            gtf_df = pd.DataFrame.from_records(
                gtf_list,
                columns=['feature', 'gene_name',
                         'low', 'high', 'strand']
            )
            try:
                gene_name = gtf_df.loc[
                    (gtf_df['feature'] == 'gene') & (
                        gtf_df['strand'] == row['strand']
                    ),
                    'gene_name'
                ].unique()[0]
                mm_gene[row['mmidx']] = gene_name
                if row['type'] == 'het_snp':
                    gene_het_snp_ratio[gene_name].append(
                        row['ratio']
                    )
            except:
                pass

    gene_allelic_ratio = defaultdict(
        lambda : default_ratio
    )
    for gene_name in gene_het_snp_ratio:
        n = len(gene_het_snp_ratio[gene_name])
        if n > 0:
            gene_allelic_ratio[gene_name] = sum(
                gene_het_snp_ratio[gene_name]
            ) / n

    mm_diff_ratio = list()
    for ri, row in mm_ratio_df.iterrows():
        mm_diff_ratio.append(
            row['ratio'] -
            gene_allelic_ratio[mm_gene[row['mmidx']]]
        )
    return mm_diff_ratio


def mm_coords_from_df(df):
    mm_coords = []
    for ri, row in df.iterrows():
        ### Reading the file
        ### Choose to split SNVs based on the dbSNP tag (last column, 0 | 1)
        chrom = row['chromosome']
        pos = row['pos']
        ref = row['ref']
        snp = row['snp']
        mm_coords.append([int(pos), base_to_number[ref], int(snp)])

    mm_coords = np.array(mm_coords)
    mm_coords = mm_coords[mm_coords[:, 0].argsort(), :]
    return mm_coords


def chrom_calculation(chrom, variables):
    ####################
    # correct strand
    read_strand = correct_read_strand(
        variables['bam_file'],
        variables['gtf_file'], chrom,
        mapq_threshold=variables['mapq_threshold'],
        gene_padding=variables['gene_padding'],
        exon_padding=variables['exon_padding']
    )

    strand = dict(
        (row['read_name'], row['read_strand'])
        for ri, row in read_strand.iterrows()
    )
    ####################
    # mismatch sites

    site_df = get_mismatch_from_cs(
        variables['bam_file'], chrom,
        mapq_threshold=variables['mapq_threshold'],
        min_allele_count=variables['min_allele_count'],
        drop_non_spliced_read=variables['drop_non_spliced_read'],
        min_dist_from_splice=variables['min_dist_from_splice']
    )
    site_df = filter_mismatch_homopoly(
        variables['genome_file'], site_df, chrom,
        homopoly_length=variables['homopoly_length'])
    site_df = filter_mismatch_simple_repeat(variables['repeat_file'], site_df, chrom)
    site_df = mark_mismatch_dbsnp(variables['snp_file'], site_df, chrom)

    if site_df.shape[0] == 0:
        return (None, None, None, None)
    ####################
    # mismatch features
    mm_coords = mm_coords_from_df(site_df)

    read_clusters = find_read_clusters(variables['bam_file'], chrom)
    mm_info_list = []
    mi_info_list = []
    for read_cluster in read_clusters:
        RC_i, RC_start, RC_end, RC_cov = read_cluster
        if RC_cov < variables['min_rc_cov']:
            continue

        mm_info = get_mm_info_from_bam(
            variables['bam_file'],
            chrom, RC_start, RC_end,
            mm_coords, strand
        )

        mm_info, ratio_dict = variant_filter(
            mm_info, chrom,
            min_AB=variables['min_AB'], min_AC=variables['min_AC'],
            min_het_snp_ratio=variables['min_het_snp_ratio'],
            max_het_snp_ratio=variables['max_het_snp_ratio'])
        mm_info_list.extend(
            ['{0}:{1}:{2}:{3}'.format(
                x, ratio_dict[x][0],
                ratio_dict[x][1], ratio_dict[x][2]
             )
             for x in mm_info]
        )
        min_n = min(variables['mi_testable_common_reads'], RC_cov)
        mi_info = mutual_information(
            mm_info, RC_i,
            mi_testable_common_reads=min_n,
            mi_testable_mono=variables['mi_testable_mono']
        )
        mi_info_list.extend(
            format_mi_output_file(mi_info)
        )

    mm_info_df = pd.DataFrame.from_records(
        [x.split(':') for x in mm_info_list],
        columns=['type',
                 'chromosome', 'pos', 'strand',
                 'change_type',
                 'read_count', 'depth', 'ratio']
    ).astype(
        {'pos': 'int', 'read_count': 'int',
         'depth': 'int', 'ratio': 'float'}
    )
    mm_info_df.loc[:, 'allelic_ratio_diff'] = diff_of_allelic_ratio(
        mm_info_df, variables['gtf_file'], chrom
    )

    mi_info_df = pd.DataFrame.from_records(
        mi_info_list,
        columns=['type',
                 'chromosome', 'pos', 'strand',
                 'change_type',
                 'mi', 'n', 'pairs', 'jakarta',
                 'mis', 'mi_cov']
    ).astype(
        {'pos': 'int', 'mi': 'float', 'n': 'int'}
    )
    return (read_strand, site_df, mm_info_df, mi_info_df)


def get_site_nearbyseq(site, genome_file):
    ntpairs = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N'}
    genome = pysam.FastaFile(genome_file)
    site_list = list()
    # +
    site_plus = site.loc[site['strand'] == '+'].copy()
    if site_plus.shape[0] > 0:
        site_plus.loc[:, 'up_seq'] = site_plus.apply(
            lambda a: genome.fetch(
                a['chromosome'], int(a['pos']) - 1, int(a['pos'])
            ),
            axis=1
        ).str.upper()
        site_plus.loc[:, 'down_seq'] = site_plus.apply(
            lambda a: genome.fetch(
                a['chromosome'], int(a['pos']) + 1, int(a['pos']) + 2
            ),
            axis=1
        ).str.upper()
        site_list.append(site_plus)
    # -
    site_minus = site.loc[site['strand'] == '-'].copy()
    if site_minus.shape[0] > 0:
        site_minus.loc[:, 'up_seq'] = site_minus.apply(
            lambda a: genome.fetch(
                a['chromosome'], int(a['pos']) + 1, int(a['pos']) + 2
            ),
            axis=1
        ).str.upper().map(ntpairs)
        site_minus.loc[:, 'down_seq'] = site_minus.apply(
            lambda a: genome.fetch(
                a['chromosome'], int(a['pos']) - 1, int(a['pos'])
            ),
            axis=1
        ).str.upper().map(ntpairs)
        site_list.append(site_minus)
    genome.close()
    return pd.concat(site_list, axis=0)


def glm_score(data, train_data):
    # data with one column label train
    train_data_dummy = pd.concat(
        [
            pd.get_dummies(
                train_data[['change_type', 'up_seq', 'down_seq']]
            ),
            train_data['allelic_ratio_diff']
        ],
        axis=1
    )
    train_label = train_data['label'].map(lambda x: 1 if x=='edit' else 0)

    lm = linear_model.LinearRegression()

    lm.fit(train_data_dummy.values, train_label.values)

    data_dummy = pd.concat(
        [
            pd.get_dummies(
                data[['change_type', 'up_seq', 'down_seq']]
            ),
            data['allelic_ratio_diff']
        ],
        axis=1
    )
    data_score = lm.predict(data_dummy.values)
    result = data.assign(
        score=data_score
    )
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="L-GIREMI (Long-read RNA-seq Genome-independent Identification of RNA Editing by Mutual Information)"
    )
    parser.add_argument(
        "-b", "--bam_file",
        help="input bam file, with cs tags, sorted and indexed",
        type=str,
        default=None,
        required=True
    )
    parser.add_argument(
        "-c", "--chromosomes",
        help = "chromosomes to be analyzed",
        nargs='*',
        type=str,
        default = ['chr1', 'chr2', 'chr3', 'chr4',
                   'chr5', 'chr6', 'chr7', 'chr8',
                   'chr9', 'chr10', 'chr11', 'chr12',
                   'chr13', 'chr14', 'chr15', 'chr16',
                   'chr17', 'chr18', 'chr19', 'chr20',
                   'chr21', 'chr22', 'chrX', 'chrY']
    )
    parser.add_argument(
        "-o", "--output_prefix",
        help = "prefix of output file",
        type=str,
        default = 'out'
    )
    parser.add_argument(
        "-t", "--thread",
        help = "cores to be used",
        type=int,
        default = 1
    )
    parser.add_argument(
        "--genome_fasta",
        help = "path of genome fasta file",
        type=str,
        default=None,
        required=True
    )
    parser.add_argument(
        "--snp_bcf",
        help = "path of dbSNP bcf file",
        type=str,
        default=None,
        required=True
    )
    parser.add_argument(
        "--repeat_txt",
        help = "path of txt file of simple repeats [chromosom, start, end] (0 based)",
        type=str,
        default=None,
        required=True
    )
    parser.add_argument(
        "--annotation_gtf",
        help = "gtf (gz and tabix indexed) file of genome annotation (gencode)",
        type=str,
        default=None,
        required=True
    )
    parser.add_argument(
        "--mapq_threshold",
        help = "Min MAPQ to be considered in bam file (default: 20)",
        type=float,
        default=20
    )
    parser.add_argument(
        "--min_allele_count",
        help = "Min allele read count (default: 2)",
        type=int,
        default=2
    )
    parser.add_argument(
        "--drop_non_spliced_read",
        help = "Drop non spliced reads (default: True)",
        type=bool,
        default=True
    )
    parser.add_argument(
        "--min_dist_from_splice",
        help = "Drop sites within the distance from splice junctions (default: 4)",
        type=int,
        default=4
    )
    parser.add_argument(
        "--gene_padding",
        help = "expand the range when searching gene gtf (default: 500)",
        type=int,
        default=500
    )
    parser.add_argument(
        "--exon_padding",
        help = "expand the range when searching exon gtf (default: 10)",
        type=int,
        default=10
    )
    parser.add_argument(
        "--min_rc_cov",
        help = "min coverage of read cluster to be considered (default: 2)",
        type=int,
        default=2
    )
    parser.add_argument(
        "--homopoly_length",
        help = "left and right sequence length to be searched for the homopoly around sites (default: 5)",
        type=int,
        default=5
    )
    parser.add_argument(
        "--min_AB",
        help = "Min mismatch ratio to be considered (default: 0.1)",
        type=float,
        default=0.1
    )
    parser.add_argument(
        "--min_AC",
        help = "Min mismatch read count to be considered (default: 3)",
        type=int,
        default=3
    )
    parser.add_argument(
        "--min_het_snp_ratio",
        help = "Min ratio to be considered as heterogenous SNPs (default: 0.35)",
        type=float,
        default=0.35
    )
    parser.add_argument(
        "--max_het_snp_ratio",
        help = "Max ratio to be considered as heterogenous SNPs (default: 0.65)",
        type=float,
        default=0.65
    )
    parser.add_argument(
        "--mi_min_common_read",
        help = "Min common read for site pairs to calculate MI (default: 6)",
        type=int,
        default=6
    )
    parser.add_argument(
        "--mi_min_read",
        help = "Min read for a variant of a site in a site pair to calculate MI (default: 1)",
        type=int,
        default=1
    )
    parser.add_argument(
        "--mip_threshold",
        help = "MI p value threshold to be used to separate RNA editing sites (default: 0.05)",
        type=float,
        default=0.05
    )
    parser.add_argument('--version', action='version', version='%(prog)s {0}'.format(__version__))
    args = parser.parse_args()

    variables = {
        'mapq_threshold': args.mapq_threshold,
        'min_allele_count': args.min_allele_count,
        'drop_non_spliced_read': args.drop_non_spliced_read,
        'min_dist_from_splice': args.min_dist_from_splice,
        'gene_padding': args.gene_padding,
        'exon_padding': args.exon_padding,
        'min_rc_cov': args.min_rc_cov,
        'homopoly_length': args.homopoly_length,
        'min_AB': args.min_AB,
        'min_AC': args.min_AC,
        'min_het_snp_ratio': args.min_het_snp_ratio,
        'max_het_snp_ratio': args.max_het_snp_ratio,
        'mi_testable_common_reads': args.mi_min_common_read,
        'mi_testable_mono': args.mi_min_read,
        'mip_threshold': args.mip_threshold,
        'bam_file': args.bam_file,
        'outprefix': args.output_prefix,
        'genome_file': args.genome_fasta,
        'snp_file': args.snp_bcf,
        'repeat_file': args.repeat_txt,
        'gtf_file': args.annotation_gtf
    }

    with mp.Pool(args.thread) as p:
        result = p.map(
            partial(chrom_calculation, variables=variables),
            args.chromosomes
        )

    # gethering data
    strand_list = list()
    site_list = list()
    mm_info_list = list()
    mi_info_list = list()
    for read_strand, site_df, mm_info_df, mi_info_df in result:
        if read_strand is not None:
            strand_list.append(read_strand)
        if site_df is not None:
            site_list.append(site_df)
        if mm_info_df is not None:
            mm_info_list.append(mm_info_df)
        if mi_info_df is not None:
            mi_info_list.append(mi_info_df)

    stranddf = pd.concat(strand_list, axis=0)
    stranddf.to_csv(
        variables['outprefix'] + '.strand',
        sep='\t', index=False
    )
    del strand_list

    sitedf = pd.concat(site_list, axis=0)
    sitedf.to_csv(
        variables['outprefix'] + '.site',
        sep='\t', index=False
    )
    del site_list

    mmdf = pd.concat(mm_info_list, axis=0)
    mmdf = get_site_nearbyseq(mmdf, variables['genome_file'])
    del mm_info_list

    # mi and mi p value
    midf = pd.concat(mi_info_list, axis=0)
    miecdf = ecdf(midf.loc[midf['type'] == 'het_snp', 'mi'])
    midf.loc[:, 'mip'] = midf['mi'].map(miecdf)
    midf.to_csv(
        variables['outprefix'] + '.mi',
        sep='\t', index=False
    )
    del mi_info_list

    # GLM score
    train_data = midf[
            ['type', 'chromosome', 'pos',
             'strand', 'change_type', 'mi', 'mip']
    ].copy()
    train_data.loc[:,'label'] = train_data.apply(
        lambda x: ['other', 'edit'][
            int(x['type'] == 'mismatch') and
            (float(x['mip']) <= variables['mip_threshold']) and
            (x['change_type'] == 'A>G')
        ],
        axis=1
    )
    train_data = pd.merge(
        train_data, mmdf,
        on=['type', 'chromosome', 'pos',
            'strand', 'change_type'],
        how='inner'
    )
    result = glm_score(mmdf, train_data)
    result.to_csv(
        variables['outprefix'] + '.score',
        sep='\t', index=False
    )