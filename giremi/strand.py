import pysam
import pandas as pd
from collections import Counter
from .cs import CS
from .utils import positions_in_intervals, merge_intervals


def get_gtf_gene_strand(gtf_gene_list, read_low, read_high, padding = 500):
    # [feature, gene_name, low, high, strand]
    if len(gtf_gene_list) == 0:
        return None
    strand_set = set([
        entry[4]
        for entry in gtf_gene_list
        if ((entry[2] - padding) <= read_high) &
        ((entry[3] + padding) >= read_low)
    ])

    if len(strand_set) == 1:
        return list(strand_set)[0]
    else:
        return None


def get_gtf_exon_strand(gtf_exon_list, read_intron_list, padding = 10):
    if len(gtf_exon_list) == 0:
        return None
    if len(read_intron_list) == 0:
        return None
    gtf_exon_list.sort(key = lambda a: a[2])
    read_intron_low = [a[0] for a in read_intron_list]
    read_intron_high = [a[1] for a in read_intron_list]
    gtf_exon_low_padding, _  = merge_intervals(
        [[a[2] - padding, a[2] + padding] for a in gtf_exon_list]
    )
    gtf_exon_high_padding, _ = merge_intervals(
        [[a[3] - padding, a[3] + padding] for a in gtf_exon_list]
    )

    read_intron_low_in_region, read_intron_low_in_region_id = \
        positions_in_intervals(read_intron_low, gtf_exon_high_padding)
    read_intron_high_in_region, read_intron_high_in_region_id = \
        positions_in_intervals(read_intron_high, gtf_exon_low_padding)
    gtf_exon_low_overlap = [
        gtf_exon_low_padding[i][0] + padding
        for i in read_intron_high_in_region_id
        if i >= 0
    ]
    gtf_exon_high_overlap = [
        gtf_exon_high_padding[i][0] + padding
        for i in read_intron_low_in_region_id
        if i >= 0
    ]

    splice_list = [
        item[4]
        for item in gtf_exon_list
        if item[2] in gtf_exon_low_overlap
    ] + [
        item[4]
        for item in gtf_exon_list
        if item[3] in gtf_exon_high_overlap
    ]

    if len(splice_list) > 0:
        read_gtf_strand_counter = list(Counter(splice_list).items())
        read_gtf_strand_counter.sort(key = lambda a: a[1], reverse = True)
        return read_gtf_strand_counter[0][0]
    else:
        return None


def get_gtf_strand(gtf_list, read_low, read_high, read_intron_list,
                   gene_padding=500, exon_padding=10):
    gtf_strand = None
    # gtf gene strand
    gtf_gene_strand = get_gtf_gene_strand(
        gtf_gene_list = [item for item in gtf_list if item[0] == 'gene'],
        read_low = read_low,
        read_high = read_high,
        padding = gene_padding
    )
    if gtf_gene_strand is not None:
        if len(gtf_gene_strand) == 1:
            gtf_strand = gtf_gene_strand[0]
        else:
            gtf_strand = get_gtf_exon_strand(
                gtf_exon_list = [item for item in gtf_list if item[0] in ['exon', 'CDS']],
                read_intron_list = read_intron_list,
                padding = exon_padding
            )
    return gtf_strand


def get_splice_strand(read_intron_list):
    splice_pattern = [
        x[3][0] + x[3][1] + x[3][-2] + x[3][-1] for x in read_intron_list
    ]
    splice_pattern_count = list(Counter(splice_pattern).items())
    splice_pattern_count.sort(key = lambda a: a[1], reverse = True)
    strand = None
    if len(splice_pattern_count) >= 1:
        pattern = splice_pattern_count[0][0]
        if pattern in ['gtag', 'gcag', 'atac']:
            strand = '+'
        elif pattern in ['ctac', 'ctgc', 'gtat']:
            strand = '-'
    return strand


def get_read_strand(seq_strand, gtf_list, read_low, read_high, read_intron_list,
                    gene_padding=500, exon_padding=10):
    # gtf strand
    gtf_strand = get_gtf_strand(
        gtf_list = gtf_list,
        read_low = read_low,
        read_high = read_high,
        read_intron_list = read_intron_list,
        gene_padding=gene_padding,
        exon_padding=exon_padding
    )
    # splice strand
    splice_strand = get_splice_strand(read_intron_list = read_intron_list)
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


def correct_read_strand_in_chrom(chrom,
                                 bam_file, gtf_file, genome_file,
                                 gene_padding = 500,
                                 exon_padding = 10,
                                 keep_non_spliced_read = False,
                                 mode = 'cs'):

    sam = pysam.AlignmentFile(bam_file, 'rb')
    gtf = pysam.TabixFile(gtf_file)
    genome = pysam.FastaFile(genome_file)

    region_gtf_list = list()
    for gtf_entry in gtf.fetch(
            chrom, parser=pysam.asGTF()
    ):
        region_gtf_list.append(
            [gtf_entry.feature,
             gtf_entry.gene_name,
             gtf_entry.start, gtf_entry.end,
             gtf_entry.strand]
        )

    read_strand_list = []

    for read in sam.fetch(chrom):
        # 0 based
        seq_strand = '-' if read.is_reverse else '+'
        if mode == 'cs':
            # mode with cs tags
            cs_tag_string = read.get_tag('cs')
            read_CS = CS.from_cs_tag_string(
                cs_tag_string,
                chrom, read.reference_start, seq_strand
            )
        else:
            # using CIGAR string and MD tags instead
            cigar_string = read.cigarstring
            md_tag_string = read.get_tag('MD')
            read_seq = read.query_sequence
            # reverse complemented for the reversed reads
            ref_seq = genome.fetch(
                chrom, read.reference_start, read.reference_end
            )
            read_CS = CS.from_cigar_string(
                cigar_string, md_tag_string, read_seq, ref_seq,
                chrom, read.reference_start, seq_strand
            )
        read_intron_list = read_CS.get_introns(coordinate = 'contig')
        read_intron_list.sort(key = lambda a: a[0])
        # ['low', 'high', 'ope', 'val']
        if not keep_non_spliced_read:
            if len(read_intron_list) == 0:
                continue
            else:
                pass
        else:
            pass
        gtf_list = [
            [g_feature, g_name, g_low, g_high, g_strand]
            for g_feature, g_name, g_low, g_high, g_strand in region_gtf_list
            if (g_low < (read.reference_start - gene_padding))
            or (g_high > (read.reference_end + gene_padding))
        ]
        # correct strand
        read_strand = get_read_strand(
            seq_strand, gtf_list = gtf_list,
            read_low = read.reference_start,
            read_high = read.reference_end,
            read_intron_list = read_intron_list,
            gene_padding = gene_padding,
            exon_padding = exon_padding
        )

        read_strand_list.append([read.query_name, read_strand])
    return chrom, read_strand_list

########################################
