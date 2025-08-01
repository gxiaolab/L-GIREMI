#! /usr/bin/env python3
import pysam
import argparse
import pandas as pd
import numpy as np
import multiprocessing as mp
from functools import partial

####################


def alu_df_aei(alu, variables):
    genome = pysam.FastaFile(variables['genome_file'])

    sam = pysam.AlignmentFile(variables['bam_file'], 'rb')

    vcf = pysam.VariantFile(variables['snp_file'])

    if variables['gtf_file'] is not None:
        gtf = pysam.TabixFile(variables['gtf_file'])
    else:
        pass

    nt_count_list = list()

    for ri, row in alu.iterrows():
        chrom = row['chromosome']
        start = int(row['start'])
        end = int(row['end'])
        name = row['name']
        alu_label = '{}:{}-{}:{}'.format(
            chrom, start, end, name
        )
        # print(alu_label)
        reads = [
            a for a in sam.fetch(chrom, start, end)
        ]
        if len(reads) == 0:
            # no reads for the region
            # print('no read')
            continue
        # print('reads: {}'.format(len(reads)))
        in_alu_snp_pos = [
            a.start for a in vcf.fetch(chrom, start, end)
        ]

        # determine main strand
        main_strand = None
        if variables['gtf_file'] is not None:
            # using gtf
            gtf_gene_strands = [
                gtf_entry.strand
                for gtf_entry in gtf.fetch(
                        chrom, start, end,
                        parser=pysam.asGTF()
                ) if gtf_entry.feature == 'gene'
            ]
            if len(gtf_gene_strands) > 0:
                strand_count = pd.Series(
                    gtf_gene_strands
                ).value_counts()
                if strand_count.shape[0] == 1:
                    main_strand = strand_count.index[0]
                elif strand_count.shape[0] > 1:
                    if strand_count['-'] != strand_count['+']:
                        main_strand = strand_count.index[
                            strand_count.argmax()
                        ]
                    else:
                        main_strand = None
                else:
                    main_strand = None

        # using read strand
        if main_strand is None:
            # read 1 main_strand
            molecule_strand1 = []
            molecule_strand2 = []
            if variables['molecule_strand'] == 'read2':
                # read1 strand is the reverse strand of the molecule
                # read2 strand is the same strand of the molecule
                molecule_strand1 = [
                    '+' if a.is_reverse else '-'
                    for a in reads if a.is_read1
                ]
                molecule_strand2 = [
                    '-' if a.is_reverse else '+'
                    for a in reads if a.is_read2
                ]
            elif variables['molecule_strand'] == 'read2':
                # read1 strand is the same strand of the molecule
                # read2 strand is the reverse strand of the molecule
                molecule_strand1 = [
                    '-' if a.is_reverse else '+'
                    for a in reads if a.is_read1
                ]
                molecule_strand2 = [
                    '+' if a.is_reverse else '-'
                    for a in reads if a.is_read2
                ]
            else:
                continue
            # print('len read_strand1: {}'.format(len(read_strand1)))
            molecule_strands = molecule_strand1 + molecule_strand2
            if len(molecule_strands) > 0:
                # print('read_strands')
                # print(read_strands)
                strand_count = pd.Series(molecule_strands).value_counts()
                # print(strand_count)
                # print(strand_count.shape)
                if strand_count.shape[0] == 1:
                    main_strand = strand_count.index[0]
                elif strand_count.shape[0] > 1:
                    if strand_count['+'] != strand_count['-']:
                        main_strand = strand_count.index[
                            strand_count.argmax()
                        ]
                    else:
                        main_strand = None
                else:
                    main_strand = None
        # print('final main_strand: {}'.format(main_strand))
        if main_strand == '+':
            nt = 'A'
        elif main_strand == '-':
            nt = 'T'
        elif main_strand is None:
            continue

        refnt = genome.fetch(chrom, start, end).upper()
        in_alu_pos = list(
            np.where([a == nt for a in list(refnt)])[0] + start
        )

        in_alu_pos_check = [
            a for a in in_alu_pos if a not in in_alu_snp_pos
        ]

        if len(in_alu_pos_check) == 0:
            continue
        dep_a, dep_c, dep_g, dep_t = sam.count_coverage(
            contig = chrom, start = start, end = end
        )

        in_alu_read_pos_check = [a - start for a in in_alu_pos_check]
        if nt == 'A':
            dep_ori = [dep_a[a] for a in in_alu_read_pos_check]
            dep_edit = [dep_g[a] for a in in_alu_read_pos_check]
        else:
            dep_ori = [dep_t[a] for a in in_alu_read_pos_check]
            dep_edit = [dep_c[a] for a in in_alu_read_pos_check]
        # calculate read and nt for each position

        if len(dep_ori) == 0:
            # no items for the alu
            continue
        region_depth = (
            sum(dep_a) + sum(dep_c) + sum(dep_g) + sum(dep_t)
        ) / len(dep_a)
        nt_count_list.append(
            [alu_label,
             region_depth,
             sum(dep_ori), sum(dep_edit)]
        )

    # calculate read aei
    if len(nt_count_list) > 0:
        aei = pd.DataFrame(
            nt_count_list,
            columns = ['alu', 'depth', 'ref', 'alt']
        )
        # calculate aei
        aei.loc[:, 'aei'] = aei['alt'] / (
            aei['alt'] + aei['ref']
        ) * 100
    else:
        aei = None
    sam.close()
    vcf.close()
    return aei

####################


def parse_args():
    parser = argparse.ArgumentParser(
        description='Calculate AEI from bam file, both in read level and Alu level'
    )
    parser.add_argument(
        '-b', '--bam_file',
        help='input bam file, with cs tags, sorted and indexed',
        type=str,
        default=None,
        required=True
    )
    parser.add_argument(
        '-c', '--chromosomes',
        help = 'chromosomes to be analyzed',
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
        '-o', '--output_prefix',
        help = 'prefix of output file',
        type=str,
        default = 'out'
    )
    parser.add_argument(
        '--genome_fasta',
        help = 'path of genome fasta file',
        type=str,
        default=None,
        required=True
    )
    parser.add_argument(
        '--snp_bcf',
        help = 'path of dbSNP bcf file',
        type=str,
        default=None,
        required=True
    )
    parser.add_argument(
        '--molecule_strand',
        help = 'molecule_strand same as read1 strand or read2 strand, read2 is default',
        type=str,
        default='read2',
        choices=['read1', 'read2', 'unstranded'],
        required=True
    )
    parser.add_argument(
        "--annotation_gtf",
        help = "gtf (gz and tabix indexed) file of genome annotation (gencode)",
        type=str,
        default=None
    )
    parser.add_argument(
        '--alu_txt',
        help = 'path of bed file of simple repeats [chromosom, start, end, name] (0 based)',
        type=str,
        default=None,
        required=True
    )
    parser.add_argument(
        '-t', '--thread',
        help = 'cores to be used',
        type=int,
        default = 1
    )
    args = parser.parse_args()
    return args


def main():
    args = parse_args()
    variables = {
        'bam_file': args.bam_file,
        'outprefix': args.output_prefix,
        'genome_file': args.genome_fasta,
        'snp_file': args.snp_bcf,
        'molecule_strand' : args.molecule_strand,
        'gtf_file' : args.annotation_gtf,
        'alu_file': args.alu_txt
    }
    alu = pd.read_table(
        variables['alu_file'], header = 0, sep = '\t'
    )
    alu = alu.loc[alu['chromosome'].isin(args.chromosomes)].copy()
    alu_df_list = [
        alu.loc[alu.index.values % args.thread == i].copy()
        for i in range(args.thread)
    ]
    # df_result = [alu_df_aei(alu, variables = variables)]
    with mp.Pool(args.thread) as p:
        df_result = p.map(
            partial(alu_df_aei, variables=variables),
            alu_df_list
        )

    aei_list = list()
    for c_aei in df_result:
        if c_aei is not None:
            aei_list.append(c_aei)
        else:
            pass
    if len(aei_list) > 0:
        aei = pd.concat(aei_list, axis = 0)
        aei.to_csv(
            variables['outprefix'] + '.aei',
            sep = '\t', index = False
        )
    else:
        pass


if __name__ == '__main__':
    main()

####################
