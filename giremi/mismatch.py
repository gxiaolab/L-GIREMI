import numpy as np
import pandas as pd
from collections import defaultdict
from .cs import CS
from .utils import merge_intervals
from .utils import positions_in_intervals
from .mutual_information import mismatch_pair_mutual_info
from .mutual_information import mean_mismatch_pair_mutual_info


def get_region_mismatches_with_filters(chromosome, start_pos, end_pos,
                                       sam, genome,
                                       keep_non_spliced_read = False,
                                       min_dist_from_splice = 4,
                                       min_allele_depth = 3,
                                       min_allele_ratio = 0.1,
                                       min_total_depth = 6,
                                       homopoly_length = 5,
                                       simple_repeat_intervals = [],
                                       snp_positions = [],
                                       read_strand_dict = None,
                                       min_het_snp_ratio = 0.35,
                                       max_het_snp_ratio = 0.65,
                                       mode = 'cs'):

    mismatches = {
        '+': defaultdict(
            lambda: {'ref': '', 'type': 'mismatch',
                     'depth': defaultdict(int),
                     'nt': defaultdict(list),
                     'up': '', 'down': ''}
        ),
        '-': defaultdict(
            lambda: {'ref': '', 'type': 'mismatch',
                     'depth': defaultdict(int),
                     'nt': defaultdict(list),
                     'up': '', 'down': ''}
        )
    }

    if read_strand_dict is None:
        read_strand_dict = dict()
    else:
        pass

    for read in sam.fetch(chromosome, start_pos, end_pos):
        # if read.mapq < mapq_threshold:
        #     continue
        # elif read.is_secondary: # skip secondary reads
        #     continue
        read_strand = '-' if read.is_reverse else '+'
        if read.query_name in read_strand_dict:
            read_strand = read_strand_dict[read.query_name]
        else:
            read_strand_dict[read.query_name] = read_strand
        if mode == 'cs':
            # mode with cs tags
            cs_tag_string = read.get_tag('cs')
            read_CS = CS.from_cs_tag_string(
                cs_tag_string,
                chromosome, read.reference_start, read_strand
            )
        else:
            # using CIGAR string and MD tags instead
            cigar_string = read.cigarstring
            md_tag_string = read.get_tag('MD')
            read_seq = read.query_sequence
            # reverse complemented for the reversed reads
            ref_seq = genome.fetch(
                chromosome, read.reference_start, read.reference_end
            )
            read_CS = CS.from_cigar_string(
                cigar_string, md_tag_string, read_seq, ref_seq,
                chromosome, read.reference_start, read_strand
            )
        read_mismatches = [
            [a[0], a[3]]
            for a in read_CS.get_mismatches(coordinate = 'contig')
        ]
        read_mismatches.sort(key = lambda a: a[0])
        # read_mismatches is a list with [ref_start_pos, "ag"]
        read_introns = read_CS.get_introns(coordinate = 'contig')
        read_introns.sort(key = lambda a: a[0])
        # read_introns is a list with
        # [ref_start_pos, ref_end_pos, '~', 'gt[length]ag']

        # check whether the are splicing in the read
        if not keep_non_spliced_read:
            if len(read_introns) == 0:
                continue
            else:
                pass
        else:
            pass
        if len(read_mismatches) > 0:
            # FILTER 1: remove mismatches close to the splicing sites.
            if len(read_introns) > 0 and min_dist_from_splice > 0:
                intron_starts = [a[0] for a in read_introns]
                intron_ends = [a[1] for a in read_introns]
                splicing_pos = intron_starts + intron_ends
                splicing_pos.sort()
                splicing_intervals, _ = merge_intervals(
                    [
                        [a - min_dist_from_splice, a + min_dist_from_splice]
                        for a in splicing_pos
                    ]
                )
                in_intervals, interval_id = positions_in_intervals(
                    [a[0] for a in read_mismatches],
                    splicing_intervals
                )
            else:
                in_intervals = [False for a in read_mismatches]

            read_mismatches_rm_splicing = [
                a for a, in_interval in zip(read_mismatches, in_intervals)
                if not in_interval
            ]
            for pos, refalt in read_mismatches_rm_splicing:
                mismatches[read_strand][pos]['ref'] = refalt[0].upper()
                mismatches[read_strand][pos]['nt'][refalt[1].upper()].append(
                    read.query_name
                )
        else:
            pass

    for strand in ['+', '-']:
        if len(mismatches[strand]) == 0:
            continue
        else:
            pass
        ####################
        # read names for ref nt
        positions = list(mismatches[strand].keys())
        pileup = sam.pileup(
            contig = chromosome, start = start_pos, stop = end_pos
        )
        for column in pileup:
            pos = column.pos
            if pos in positions:
                ref = mismatches[strand][pos]['ref']
                raw_read_names = column.get_query_names()
                raw_read_seqs = [
                    a.upper() for a in column.get_query_sequences()
                ]
                seq_keep_idx = [
                    True if a == ref else False
                    for a in raw_read_seqs
                ]
                # should check read strand
                strand_keep_idx = [
                    True if read_strand_dict[rn] == strand else False
                    for rn in raw_read_names
                ]
                read_names = [
                    a
                    for a, keep1, keep2
                    in zip(raw_read_names, seq_keep_idx, strand_keep_idx)
                    if keep1 and keep2
                ]
                mismatches[strand][pos]['nt'][ref].extend(read_names)
            else:
                pass
        ####################
        # depth for mismatches
        positions = list(mismatches[strand].keys())
        # region depth
        # r_a, r_c, r_g, r_t = sam.count_coverage(
        #     contig = chromosome, start = start_pos, stop = end_pos
        # )
        # for pos in positions:
        #     mismatches[pos]['depth']['a'] = r_a[pos - start_pos]
        #     mismatches[pos]['depth']['c'] = r_c[pos - start_pos]
        #     mismatches[pos]['depth']['g'] = r_g[pos - start_pos]
        #     mismatches[pos]['depth']['t'] = r_t[pos - start_pos]
        for pos in positions:
            nt_seqs = list(mismatches[strand][pos]['nt'].keys())
            for nt in nt_seqs:
                mismatches[strand][pos]['depth'][nt] = len(
                    mismatches[strand][pos]['nt'][nt]
                )
        ####################
        # Filters
        # FILTER 2: remove alleles of mismatches with too shallow depth
        positions = list(mismatches[strand].keys())
        for pos in positions:
            nt_seqs = list(mismatches[strand][pos]['nt'].keys())
            for nt in nt_seqs:
                a_depth = mismatches[strand][pos]['depth'][nt]
                if a_depth < min_allele_depth:
                    mismatches[strand][pos]['nt'].pop(nt)
                else:
                    pass
        # FILTER 3: remove alleles of mismatches with too small allelic ratio
        positions = list(mismatches[strand].keys())
        for pos in positions:
            t_depth = sum(mismatches[strand][pos]['depth'].values())
            nt_seqs = list(mismatches[strand][pos]['nt'].keys())
            for nt in nt_seqs:
                a_depth = mismatches[strand][pos]['depth'][nt]
                a_ratio = a_depth / t_depth
                if a_ratio < min_allele_ratio:
                    mismatches[strand][pos]['nt'].pop(nt)
                else: pass
        # FILTER 4: remove mismatches with too shallow total depth
        positions = list(mismatches[strand].keys())
        for pos in positions:
            t_depth = sum(mismatches[strand][pos]['depth'].values())
            if t_depth < min_total_depth:
                mismatches[strand].pop(pos)
            else:
                pass
        # FILTER 5: remove mismatches without enought allele after previous fitlers
        positions = list(mismatches[strand].keys())
        for pos in positions:
            nt_count = len(mismatches[strand][pos]['nt'])
            if nt_count < 2:
                mismatches[strand].pop(pos)
            else:
                pass
        # FILTER 6: remove mismatches in homopoly and assign up and down seq
        positions = list(mismatches[strand].keys())
        for pos in positions:
            left = genome.fetch(
                chromosome, pos - homopoly_length, pos
            ).upper()
            mismatches[strand][pos]['up'] = left[-1].upper()
            right = genome.fetch(
                chromosome, pos + 1, pos + homopoly_length + 1
            ).upper()
            mismatches[strand][pos]['down'] = right[0].upper()
            in_homo = (len(set(left)) == 1) or \
                (len(set(right)) == 1) or \
                (len(set(left[-int(homopoly_length / 2):] + right[0:int(homopoly_length / 2)])) == 1)
            if in_homo:
                mismatches[strand].pop(pos)
            else:
                pass
        # FILTER 7: remove mismatches in simple repeats
        positions = list(mismatches[strand].keys())
        in_repeats, interval_id = positions_in_intervals(
            positions,
            simple_repeat_intervals
        )
        positions_to_remove = [a for a, in_repeat in zip(positions, in_repeats) if in_repeat]
        for pos in positions_to_remove:
            mismatches[strand].pop(pos)
        ####################
        # update depth for mismatches
        positions = list(mismatches[strand].keys())
        for pos in positions:
            nt_seqs = list(mismatches[strand][pos]['depth'].keys())
            for nt in nt_seqs:
                mismatches[strand][pos]['depth'].pop(nt)
            nt_seqs = list(mismatches[strand][pos]['nt'].keys())
            for nt in nt_seqs:
                mismatches[strand][pos]['depth'][nt] = len(mismatches[strand][pos]['nt'][nt])
        ####################
        # Mark SNPs: mark the type of the mismatches as het_snp, snp, or mismatch
        positions = list(mismatches[strand].keys())
        for pos in positions:
            if pos in snp_positions:
                t_depth = sum(mismatches[strand][pos]['depth'].values())
                a_depth = mismatches[strand][pos]['depth'].values()
                a_ratio = [a / t_depth for a in a_depth]
                a_ratio.sort()
                major_allele_ratio = a_ratio[-1]
                # minor_allele_ratio = a_ratio[-2]
                if major_allele_ratio >= min_het_snp_ratio and major_allele_ratio <= max_het_snp_ratio:
                    mismatches[strand][pos]['type'] = 'het_snp'
                else:
                    mismatches[strand][pos]['type'] = 'snp'
            else:
                pass
    ####################
    return mismatches


def region_mismatch_analysis(chromosome, start_pos, end_pos,
                             sam, genome,
                             keep_non_spliced_read = False,
                             min_dist_from_splice = 4,
                             min_allele_depth = 3,
                             min_allele_ratio = 0.1,
                             min_total_depth = 6,
                             homopoly_length = 5,
                             simple_repeat_intervals = [],
                             snp_positions = [],
                             read_strand_dict = None,
                             min_het_snp_ratio = 0.35,
                             max_het_snp_ratio = 0.65,
                             mode = 'cs',
                             min_common_reads = 5):

    mismatches = get_region_mismatches_with_filters(
        chromosome = chromosome, start_pos = start_pos, end_pos = end_pos,
        sam = sam, genome = genome,
        keep_non_spliced_read = keep_non_spliced_read,
        min_dist_from_splice = min_dist_from_splice,
        min_allele_depth = min_allele_depth,
        min_allele_ratio = min_allele_ratio,
        min_total_depth = min_total_depth,
        homopoly_length = homopoly_length,
        simple_repeat_intervals = simple_repeat_intervals,
        snp_positions = snp_positions,
        read_strand_dict = read_strand_dict,
        min_het_snp_ratio = min_het_snp_ratio,
        max_het_snp_ratio = max_het_snp_ratio,
        mode = mode
    )

    mismatch_pair_mi_full = {'+': [], '-': []}
    mismatch_pair_mi = {'+': [], '-': []}
    mismatch_mean_mi = {'+': [], '-': []}
    for strand in ['+', '-']:
        if len(mismatches[strand]) > 1:
            mismatch_pair_mi_full[strand] = mismatch_pair_mutual_info(
                mismatches[strand], min_common_reads = min_common_reads
            )
            if len(mismatch_pair_mi_full[strand]) > 0:
                mismatch_pair_mi[strand] = [
                    a for a in mismatch_pair_mi_full[strand]
                    if a[1] == 'het_snp' or a[3] == 'het_snp'
                ]
                if len(mismatch_pair_mi[strand]) > 0:
                    mismatch_mean_mi[strand] = mean_mismatch_pair_mutual_info(
                        mismatch_pair_mi[strand]
                    )
            else:
                pass
        else:
            pass

    # generate pands DataFrame for output
    df_mismatch_pair_mi = pd.DataFrame.from_records(
        [
            [chromosome, '+', p1, type1, p2, type2, p1p2mi]
            for p1, type1, p2, type2, p1p2mi in mismatch_pair_mi['+']
        ] + [
            [chromosome, '-', p1, type1, p2, type2, p1p2mi]
            for p1, type1, p2, type2, p1p2mi in mismatch_pair_mi['-']
        ],
        columns = ['chromosome', 'strand',
                   'site1_pos', 'site1_type',
                   'site2_pos', 'site2_type', 'mi']
    )
    mismatch_records = []
    nt_pair = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    for strand in ['+', '-']:
        mean_mi_dict = dict(mismatch_mean_mi[strand])
        # calculate regional allelic ratio
        het_snp_ratios = list()
        for pos in mismatches[strand]:
            postype = mismatches[strand][pos]['type']
            if postype == 'het_snp':
                depth = mismatches[strand][pos]['depth'].copy()
                ref_allele = mismatches[strand][pos]['ref']
                total_depth = sum(depth[nt] for nt in depth)
                # ref_depth = mismatches[strand][pos]['depth'][ref_allele]
                alt_allele_depth = [
                    [nt, depth[nt]] for nt in depth if nt != ref_allele
                ]
                alt_allele_depth.sort(key = lambda a: a[1], reverse = True)
                alt_major_allele = alt_allele_depth[0][0]
                alt_major_allele_depth = alt_allele_depth[0][1]
                alt_major_ratio = alt_major_allele_depth / total_depth
                het_snp_ratios.append(alt_major_ratio)
            else:
                pass
        if len(het_snp_ratios) > 0:
            allelic_ratio = sum(het_snp_ratios) / len(het_snp_ratios)
        else:
            allelic_ratio = 0.5
        for pos in mismatches[strand]:
            depth = mismatches[strand][pos]['depth'].copy()
            postype = mismatches[strand][pos]['type']
            ref_allele = mismatches[strand][pos]['ref']
            total_depth = sum(depth[nt] for nt in depth)
            # ref_depth = mismatches[strand][pos]['depth'][ref_allele]
            alt_allele_depth = [
                [nt, depth[nt]] for nt in depth if nt != ref_allele
            ]
            alt_allele_depth.sort(key = lambda a: a[1], reverse = True)
            alt_major_allele = alt_allele_depth[0][0]
            alt_major_allele_depth = alt_allele_depth[0][1]
            alt_major_ratio = alt_major_allele_depth / total_depth
            allelic_ratio_diff = alt_major_ratio - allelic_ratio
            nt_depth = '{}:{}:{}:{}'.format(depth['A'], depth['C'], depth['T'], depth['G'])
            change_type = ''
            if strand == '+':
                change_type = '{}>{}'.format(ref_allele, alt_major_allele)
            else:
                change_type = '{}>{}'.format(
                    nt_pair[ref_allele], nt_pair[alt_major_allele]
                )
            up_nt = ''
            down_nt = ''
            if strand == '+':
                up_nt = mismatches[strand][pos]['up']
                down_nt = mismatches[strand][pos]['down']
            else:
                up_nt = nt_pair[mismatches[strand][pos]['down']]
                down_nt = nt_pair[mismatches[strand][pos]['up']]
            if pos in mean_mi_dict:
                mean_mi = mean_mi_dict[pos]
            else:
                mean_mi = np.nan
            mismatch_records.append(
                [postype, chromosome, strand, pos,
                 ref_allele, change_type, alt_major_ratio, allelic_ratio_diff,
                 total_depth, nt_depth, up_nt, down_nt, mean_mi]
            )

    df_mismatches = pd.DataFrame.from_records(
        mismatch_records,
        columns = ['type', 'chromosome', 'strand', 'pos', 'ref', 'change_type',
                   'ratio', 'allelic_ratio_diff', 'depth', 'A:C:T:G',
                   'up_seq', 'down_seq', 'mean_mi']
    )
    return df_mismatches, df_mismatch_pair_mi


########################################
