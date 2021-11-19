#! /usr/bin/env python3
import pysam
import argparse
import pandas as pd
import numpy as np
import multiprocessing as mp
from functools import partial

####################


def chrom_aei(chrom, variables):
    genome = pysam.FastaFile(variables['genome_file'])

    alu = pd.read_table(
        variables['alu_file'], header = 0, sep = '\t'
    )
    alu = alu.loc[alu['chromosome'] == chrom]

    strand = pd.read_table(variables['strand_file'], header = 0, sep = '\t')
    strand.index = strand['read_name']

    sam = pysam.AlignmentFile(variables['bam_file'], 'rb')

    vcf = pysam.VariantFile(variables['snp_file'])

    read_count_df_list = list()

    for ri, row in alu.iterrows():
        chrom = row['chromosome']
        start = row['start']
        end = row['end']
        name = row['name']
        alu_label = '{}:{}-{}:{}'.format(
            chrom, start, end, name
        )
        in_alu_snp_list = [a for a in vcf.fetch(chrom, start, end)]

        in_alu_snp_pos = [a.start for a in in_alu_snp_list]

        read_name = [
            a.qname for a in sam.fetch(chrom, start, end)
            if a.qname in strand['read_name']
        ]

        if len(read_name) == 0:
            continue

        main_strand = strand.loc[
            read_name, 'read_strand'
        ].value_counts().sort_values(ascending=False).index[0]

        if main_strand == '+':
            nt = 'A'
            nt_edit = 'G'
        elif main_strand == '-':
            nt = 'T'
            nt_edit = 'C'

        refnt = genome.fetch(chrom, start, end).upper()
        in_alu_pos = list(
            np.where([a == nt for a in list(refnt)])[0] + start
        )

        in_alu_pos_check = [
            a for a in in_alu_pos if a not in in_alu_snp_pos
        ]

        if len(in_alu_pos_check) == 0:
            continue
        pos_df_list = list()
        # calculate read and nt for each position
        for pileupcolumn in sam.pileup(contig = chrom,
                                       start = start,
                                       end = end):
            pos = pileupcolumn.pos
            if pos in in_alu_pos_check:
                raw_read_names = pileupcolumn.get_query_names()
                raw_read_seqs = [
                    a.upper() for a in pileupcolumn.get_query_sequences()
                ]
                keep_idx = [
                    True if a in [nt, nt_edit] else False
                    for a in raw_read_seqs
                ]
                read_names = [
                    a for a, b in zip(raw_read_names, keep_idx) if b
                ]
                read_seqs = [
                    a for a, b in zip(raw_read_seqs, keep_idx) if b
                ]
                if len(read_names) > 0:
                    pos_df = pd.DataFrame(
                        {
                            'alu':[alu_label] * len(read_names),
                            'pos': [pos] * len(read_names),
                            'read_name': read_names,
                            'nt': [
                                'ref' if a == nt else 'alt'
                                for a in read_seqs
                            ]
                        }
                    )
                    pos_df_list.append(pos_df)

        if len(pos_df_list) == 0:
            # no items for the alu
            continue

        alu_df = pd.concat(pos_df_list)
        if alu_df.shape[0] == 0:
            # no items for the alu, redundant
            continue
        read_count_df = alu_df.groupby(
            ['alu', 'read_name']
        )['nt'].value_counts()

        read_count_df.name = 'count'

        read_count_df_list.append(
            read_count_df.reset_index(drop = False)
        )

    # calculate read aei
    try:
        read_aei = pd.concat(
            read_count_df_list
        ).set_index(
            ['alu', 'read_name', 'nt']
        ).unstack(fill_value=0).reset_index()
        # rename the columns
        colname = [read_aei.columns[0][0],
                   read_aei.columns[1][0]]
        for i in range(2, len(read_aei.columns)):
            colname += [read_aei.columns[i][1]]
        read_aei.columns = colname
        for x in ['ref', 'alt']:
            if x not in read_aei.columns:
                read_aei.loc[:, x] = 0
        # calculate AEI
        read_aei.loc[:, 'aei'] = read_aei['alt'] / (
            read_aei['alt'] + read_aei['ref']
        ) * 100

        # calculate aei
        aei = read_aei[
            ['alu', 'read_name', 'alt', 'ref']
        ].groupby('alu')[['alt', 'ref']].sum().reset_index()

        aei.loc[:, 'aei'] = aei['alt'] / (
            aei['alt'] + aei['ref']
        ) * 100
    except ValueError:
        print(chrom + ' has no data output.')
        read_aei = None
        aei = None
    sam.close()
    vcf.close()
    return read_aei, aei

####################


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Calculate AEI from bam file, both in read level and Alu level"
    )
    parser.add_argument(
        "-b", "--bam_file",
        help="input bam file, with cs tags, sorted and indexed",
        type=str,
        default=None,
        required=True
    )
    parser.add_argument(
        "--strand_file",
        help = "corrrect strand file for reads, generated by the L-GIREMI",
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
        "--alu_txt",
        help = "path of bed file of simple repeats [chromosom, start, end, name] (0 based)",
        type=str,
        default=None,
        required=True
    )
    parser.add_argument(
        "-t", "--thread",
        help = "cores to be used",
        type=int,
        default = 1
    )
    args = parser.parse_args()

    variables = {
        'bam_file': args.bam_file,
        'strand_file': args.strand_file,
        'outprefix': args.output_prefix,
        'genome_file': args.genome_fasta,
        'snp_file': args.snp_bcf,
        'alu_file': args.alu_txt
    }

    with mp.Pool(args.thread) as p:
        chrom_result = p.map(
            partial(chrom_aei, variables=variables),
            args.chromosomes
        )

    read_aei_list = list()
    aei_list = list()
    for c_read_aei, c_aei in chrom_result:
        if c_read_aei is not None:
            read_aei_list.append(c_read_aei)
        if c_aei is not None:
            aei_list.append(c_aei)

    read_aei = pd.concat(read_aei_list, axis = 0)
    read_aei.to_csv(
        variables['outprefix'] + '.read_aei',
        sep = '\t', index = False
    )
    aei = pd.concat(aei_list, axis = 0)
    aei.to_csv(
        variables['outprefix'] + '.aei',
        sep = '\t', index = False
    )

####################
