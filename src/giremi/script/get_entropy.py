#! /usr/bin/env python3
import pysam
import argparse
import pandas as pd
import numpy as np
import multiprocessing as mp
from functools import partial

####################


def chrom_entropy(chrom, variables):
    genome = pysam.FastaFile(variables['genome_file'])

    region = pd.read_table(
        variables['region_file'], header = 0, sep = '\t'
    )
    region = region.loc[region['chromosome'] == chrom]

    strand = pd.read_table(variables['strand_file'], header = 0, sep = '\t')
    strand.index = strand['read_name']

    sam = pysam.AlignmentFile(variables['bam_file'], 'rb')

    vcf = pysam.VariantFile(variables['snp_file'])

    read_count_df_list = list()

    for ri, row in region.iterrows():
        chrom = row['chromosome']
        start = row['start']
        end = row['end']
        name = row['name']
        region_label = '{}:{}-{}:{}'.format(
            chrom, start, end, name
        )

        in_region_snp_pos = [
            a.start for a in vcf.fetch(chrom, start, end)
        ]

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
        in_region_pos_all = list(
            np.where([a == nt for a in list(refnt)])[0] + start
        )

        in_region_pos = [
            a for a in in_region_pos_all
            if a not in in_region_snp_pos
        ]

        if len(in_region_pos) == 0:
            continue
        pos_df_list = list()
        # calculate read and nt for each position
        for pileupcolumn in sam.pileup(contig = chrom,
                                       start = start,
                                       end = end):
            pos = pileupcolumn.pos
            if pos in in_region_pos:
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
                            'region': [region_label] * len(read_names),
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

        region_df = pd.concat(pos_df_list)
        if region_df.shape[0] == 0:
            # no items for the alu, redundant
            continue
        read_count_df = region_df.groupby(
            ['region', 'read_name']
        )['nt'].value_counts()

        read_count_df.name = 'count'

        read_count_df_list.append(
            read_count_df.reset_index(drop = False)
        )

    # calculate read statistic
    if len(read_count_df_list) > 0:
        read_nt = pd.concat(
            read_count_df_list
        ).set_index(
            ['region', 'read_name', 'nt']
        ).unstack(fill_value=0).reset_index()
        # rename the columns
        colname = [read_nt.columns[0][0],
                   read_nt.columns[1][0]]
        for i in range(2, len(read_nt.columns)):
            colname += [read_nt.columns[i][1]]
        read_nt.columns = colname
        for x in ['ref', 'alt']:
            if x not in read_nt.columns:
                read_nt.loc[:, x] = 0
            else:
                pass
        # calculate ratio
        read_nt.loc[:, 'ratio'] = read_nt['alt'] / (
            read_nt['alt'] + read_nt['ref']
        )
        read_nt.loc[:, 'ratio_bin'] = pd.cut(
            read_nt['ratio'],
            bins = [-0.1, 0.2, 0.4, 0.6, 0.8, 1]
        )
        # calculate entropy
        region_count = read_nt[
            ['region', 'read_name']
        ].groupby(
            ['region']
        )['read_name'].count().reset_index()
        ratio_bin = pd.merge(
            read_nt[
                ['region', 'read_name', 'ratio_bin']
            ].groupby(
                ['region', 'ratio_bin']
            )['read_name'].count(
            ).reset_index(),
            region_count,
            how = 'inner',
            left_on = 'region',
            right_on = 'region'
        )
        ratio_bin.columns = [
            'region', 'ratio_bin',
            'bin_rc', 'region_rc'
        ]
        ratio_bin = ratio_bin.loc[
            ratio_bin['bin_rc'] > 0
        ]
        ratio_bin.loc[:, 'p'] = ratio_bin['bin_rc'] / ratio_bin['region_rc']
        ratio_bin.loc[:, 'logp'] = np.log10(
            ratio_bin.loc[:, 'p']
        )
        ratio_bin.loc[:, 'plogp'] = ratio_bin['p'] * ratio_bin['logp']
        entropy = ratio_bin[
            ['region', 'region_rc', 'plogp']
        ].groupby(
            ['region', 'region_rc']
        )['plogp'].sum().reset_index()

        entropy.loc[:, 'entropy'] = - entropy['plogp']
        del entropy['plogp']
        entropy.columns = ['region', 'region_read_count', 'entropy']
    else:
        entropy = None
    sam.close()
    vcf.close()
    return entropy

####################


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Calculate Shannon entropy for regions from bam file"
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
        "--region_txt",
        help = "path of bed file of regions [chromosom, start, end, name] (0 based)",
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
        'region_file': args.region_txt
    }

    with mp.Pool(args.thread) as p:
        chrom_result = p.map(
            partial(chrom_entropy,
                    variables=variables),
            args.chromosomes
        )

    entropy_list = list()
    for c_entropy in chrom_result:
        if c_entropy is not None:
            entropy_list.append(c_entropy)
        else:
            pass
    if len(entropy_list) > 0:
        entropy = pd.concat(
            entropy_list, axis = 0
        )
        entropy.to_csv(
            variables['outprefix'] + '.entropy',
            sep = '\t', index = False
        )
    else:
        pass

####################
