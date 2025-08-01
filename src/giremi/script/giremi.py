#! /usr/bin/env python3
import argparse
import pysam
import logging
import pandas as pd
import numpy as np
import multiprocessing as mp
from functools import partial
import giremi
from giremi.footprint import get_footprints
from giremi.fileio import read_simple_repeat_intervals
from giremi.fileio import read_snp_positions_in_region
from giremi.mismatch import region_mismatch_analysis
from giremi.strand import correct_read_strand_in_region
from giremi.stat import ecdf
from giremi.stat import glm_score
from giremi.stat import score_performance


def footprint_bulk_calculation(footprints, variables):
    sam = pysam.AlignmentFile(variables['bam_file'], 'rb')
    genome = pysam.FastaFile(variables['genome_file'])
    vcf = pysam.VariantFile(variables['snp_file'])
    gtf = pysam.TabixFile(variables['gtf_file'])

    repeats = read_simple_repeat_intervals(variables['repeat_file'])

    mismatch_df_list = []
    mi_df_list = []
    strand_list = []
    removed_mismatch_df_list = []
    for footprint in footprints:
        chromosome, start_pos, end_pos, rc = footprint
        # correct read strand
        if not variables['skip_strand_correction']:
            read_strand_list = correct_read_strand_in_region(
                chromosome, start_pos, end_pos,
                sam, gtf, genome,
                variables['gene_padding'], variables['exon_padding'],
                keep_non_spliced_read = variables['keep_non_spliced_read'],
                mode = variables['mode']
            )
            strand_list.extend(read_strand_list)
            read_strand_dict = dict([
                [rname, correctedstrand]
                for rname, oldstrand, correctedstrand in read_strand_list
            ])
        else:
            read_strand_dict = None
        # get snps in the footprint
        snp_positions = read_snp_positions_in_region(
            vcf, chromosome, start_pos, end_pos
        )
        # get simple repeats in the footprint
        simple_repeat_intervals = [
            [rs, re]
            for rs, re in repeats[chromosome]
            if (rs > end_pos) or (re < start_pos)
        ]
        df_mismatches, df_mis_pair_mi, df_removed_mismatches = region_mismatch_analysis(
            chromosome, start_pos, end_pos, sam, genome,
            simple_repeat_intervals=simple_repeat_intervals,
            snp_positions=snp_positions,
            keep_non_spliced_read=variables['keep_non_spliced_read'],
            min_dist_from_splice=variables['min_dist_from_splice'],
            min_allele_depth=variables['min_allele_depth'],
            min_allele_ratio=variables['min_allele_ratio'],
            min_total_depth=variables['min_total_depth'],
            homopoly_length=variables['homopoly_length'],
            min_het_snp_ratio=variables['min_het_snp_ratio'],
            max_het_snp_ratio=variables['max_het_snp_ratio'],
            mismatch_window_size=variables['mismatch_window_size'],
            max_window_mismatch=variables['max_window_mismatch'],
            max_window_mismatch_type=variables['max_window_mismatch_type'],
            min_common_reads=variables['mi_min_common_reads'],
            mode=variables['mode'],
            read_strand_dict=read_strand_dict
        )
        mismatch_df_list.append(df_mismatches)
        mi_df_list.append(df_mis_pair_mi)
        removed_mismatch_df_list.append(df_removed_mismatches)
    mismatch_df = pd.concat(mismatch_df_list, axis=0)
    mi_df = pd.concat(mi_df_list, axis=0)
    strand_df = pd.DataFrame.from_records(
        strand_list,
        columns = ['read_name', 'original_read_strand', 'corrected_read_strand']
    )
    removed_mismatch_df = pd.concat(removed_mismatch_df_list, axis=0)
    sam.close()
    genome.close()
    vcf.close()
    gtf.close()
    return mismatch_df, mi_df, strand_df, removed_mismatch_df


def run_glm(mmdf, mip_threshold=0.05, model='lm'):
    train_data_pos = mmdf.loc[
        mmdf['mean_mi'].notna() &
        (mmdf['mip'] <= mip_threshold) &
        (mmdf['type'] == 'mismatch'),
        ['type', 'chromosome', 'pos',
         'strand', 'change_type', 'mean_mi', 'mip']
    ].copy().reset_index(drop=True)
    train_data_pos.loc[:, 'label'] = 'train_edit'
    pos_num = train_data_pos.shape[0]
    train_data_neg = mmdf.loc[
        mmdf['mean_mi'].notna() &
        (mmdf['mip'] > mip_threshold) &
        (mmdf['type'] != 'mismatch'),
        ['type', 'chromosome', 'pos',
         'strand', 'change_type', 'mean_mi', 'mip']
    ].copy().reset_index(drop=True)
    neg_num = train_data_neg.shape[0]
    train_data_neg.loc[:, 'label'] = 'train_other'
    bnum = min(pos_num, neg_num)
    train_data = pd.concat([
        train_data_pos.sample(n=bnum, replace=False),
        train_data_neg.sample(n=bnum, replace=False)
    ])
    mmdf_labeled = pd.merge(
        mmdf,
        train_data[
            ['type', 'chromosome', 'pos', 'strand', 'change_type', 'label']
        ],
        on = ['type', 'chromosome', 'pos',
              'strand', 'change_type'],
        how = 'left'
    )
    # select functions according to features
    glmresult = glm_score(mmdf_labeled, model)

    scorepf = score_performance(
        tp = train_data_pos.copy(),
        tn = train_data_neg.copy(),
        score = glmresult
    )
    return glmresult, scorepf


def parse_args():
    parser = argparse.ArgumentParser(
        description='L-GIREMI (Long-read RNA-seq Genome-independent Identification of RNA Editing by Mutual Information)'
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
        help='chromosomes to be analyzed',
        nargs='*',
        type=str,
        default=['chr1', 'chr2', 'chr3', 'chr4',
                 'chr5', 'chr6', 'chr7', 'chr8',
                 'chr9', 'chr10', 'chr11', 'chr12',
                 'chr13', 'chr14', 'chr15', 'chr16',
                 'chr17', 'chr18', 'chr19', 'chr20',
                 'chr21', 'chr22', 'chrX', 'chrY']
    )
    parser.add_argument(
        '-o', '--output_prefix',
        help='prefix of output file',
        type=str,
        default='out'
    )
    parser.add_argument(
        '-t', '--thread',
        help='cores to be used',
        type=int,
        default = 1
    )
    parser.add_argument(
        '--annotation_gtf',
        help='gtf (gz and tabix indexed) file of genome annotation (gencode)',
        type=str,
        default=None,
        required=True
    )
    parser.add_argument(
        '--genome_fasta',
        help='path of genome fasta file',
        type=str,
        default=None,
        required=True
    )
    parser.add_argument(
        '--homopoly_length',
        help='left and right sequence length to be searched for the homopoly around sites (default: 5)',
        type=int,
        default=5
    )
    parser.add_argument(
        '--keep_non_spliced_read',
        help='Keep non spliced reads (default: without this parameter, reads without splicing would be dropped.)',
        action='store_true'
    )
    parser.add_argument(
        '--max_het_snp_ratio',
        help='Max ratio to be considered as heterogenous SNPs (default: 0.65)',
        type=float,
        default=0.65
    )
    parser.add_argument(
        '--mi_calculation_only',
        help='Only calculate the mutual information, without the modeling step. (default: calculate the steps after the MI step)',
        action='store_true'
    )
    parser.add_argument(
        '--mi_min_common_read',
        help='Min common read for site pairs to calculate MI (default: 6)',
        type=int,
        default=6
    )
    parser.add_argument(
        '--mi_p_threshold',
        help='MI p value threshold to be used to separate RNA editing sites (default: 0.05)',
        type=float,
        default=0.05
    )
    parser.add_argument(
        '--min_allele_depth',
        help='Min mismatch allele read count to be considered (default: 3)',
        type=int,
        default=3
    )
    parser.add_argument(
        '--min_allele_ratio',
        help='Min mismatch ratio to be considered (default: 0.05)',
        type=float,
        default=0.05
    )
    parser.add_argument(
        '--min_het_snp_ratio',
        help='Min ratio to be considered as heterogenous SNPs (default: 0.35)',
        type=float,
        default=0.35
    )
    parser.add_argument(
        '--min_total_depth',
        help='Min mismatch total read count to be considered (default: 6)',
        type=int,
        default=2
    )
    parser.add_argument(
        '--min_dist_from_splice',
        help='Drop sites within the distance from splice junctions (default: 4)',
        type=int,
        default=4
    )
    parser.add_argument(
        '--mismatch_window_size',
        help='Window length to search the number of mismatches, too many mismatches indicate a region of many sequencing errors.',
        type=int,
        default=100
    )
    parser.add_argument(
        '--max_window_mismatch',
        help='Too many mismaches indicates a region of many sequencing errors.',
        type=int,
        default=10
    )
    parser.add_argument(
        '--max_window_mismatch_type',
        help='Too many mismaches indicates a region of many sequencing errors.',
        type=int,
        default=3
    )
    parser.add_argument(
        '--mode',
        help='Mode to find mismatches in BAM: with cs tags (cs), or CIGAR and MD tags (cigar). (defalt: cs)',
        type=str,
        default = 'cs'
    )
    parser.add_argument(
        '--model',
        help='GLM model for scoring: lm, logistic. (defalt: lm)',
        type=str,
        default = 'lm'
    )
    parser.add_argument(
        '--padding_exon',
        help='expand the range when searching exon gtf (default: 10)',
        type=int,
        default=10
    )
    parser.add_argument(
        '--padding_gene',
        help='expand the range when searching gene gtf (default: 500)',
        type=int,
        default=500
    )
    parser.add_argument(
        '--repeat_txt',
        help='path of txt file of simple repeats [chromosom, start, end] (0 based)',
        type=str,
        default=None,
        required=True
    )
    parser.add_argument(
        '--skip_strand_correction',
        help='skip the strand correction step to save runtime.',
        action='store_true'
    )
    parser.add_argument(
        '--snp_bcf',
        help='path of dbSNP bcf file',
        type=str,
        default=None,
        required=True
    )
    VERSION = '0.2.4'
    parser.add_argument(
        '--version',
        action='version',
        version='%(prog)s {0}'.format(giremi.__version__)
    )
    args = parser.parse_args()
    return args


def main():
    args = parse_args()

    variables = {
        'bam_file' : args.bam_file,
        'genome_file' : args.genome_fasta,
        'gtf_file' : args.annotation_gtf,
        'repeat_file' : args.repeat_txt,
        'snp_file' : args.snp_bcf,
        'exon_padding' : args.padding_exon,
        'gene_padding' : args.padding_gene,
        'homopoly_length' : args.homopoly_length,
        'keep_non_spliced_read' : args.keep_non_spliced_read,
        'max_het_snp_ratio' : args.max_het_snp_ratio,
        'min_allele_depth' : args.min_allele_depth,
        'min_allele_ratio' : args.min_allele_ratio,
        'min_dist_from_splice' : args.min_dist_from_splice,
        'min_het_snp_ratio' : args.min_het_snp_ratio,
        'min_total_depth' : args.min_total_depth,
        'mip_threshold' : args.mi_p_threshold,
        'mi_min_common_reads' : args.mi_min_common_read,
        'mi_calculation_only': args.mi_calculation_only,
        'mismatch_window_size': args.mismatch_window_size,
        'max_window_mismatch': args.max_window_mismatch,
        'max_window_mismatch_type': args.max_window_mismatch_type,
        'mode': args.mode,
        'model': args.model,
        'skip_strand_correction': args.skip_strand_correction,
    }

    # logging

    logging.basicConfig(
        format='%(asctime)s %(levelname)s %(message)s',
        level=logging.INFO
    )

    # footprint depth use variables['min_total_depth']
    message = 'Get regions that are covered by enough reads.'
    logging.info(message)
    footprints = get_footprints(
        variables['bam_file'], args.chromosomes, variables['min_total_depth']
    )
    n = int(len(footprints) / args.thread / 2)
    footprint_chunk = []
    for i in range(0, len(footprints), n):
        footprint_chunk.append(footprints[i:(i + n)])

    # calculate mismaches
    message = 'Calculate mismatches in each region.'
    logging.info(message)
    with mp.Pool(args.thread) as p:
        footprint_result = p.map(
            partial(footprint_bulk_calculation,
                    variables=variables),
            footprint_chunk
        )
    mismatch_df_list = []
    mi_df_list = []
    strand_df_list = []
    removed_mismatch_df_list = []
    for dfmismatch, dfmi, dfstrand, dfremoved in footprint_result:
        mismatch_df_list.append(dfmismatch)
        mi_df_list.append(dfmi)
        strand_df_list.append(dfstrand)
        removed_mismatch_df_list.append(dfremoved)

    mismatch_df = pd.concat(mismatch_df_list, axis = 0)
    mi_df = pd.concat(mi_df_list, axis = 0)
    strand_df = pd.concat(strand_df_list, axis = 0)
    removed_mismatch_df = pd.concat(removed_mismatch_df_list, axis = 0)
    # output strand table
    strand_df.to_csv(
        args.output_prefix + '.strand.txt',
        sep = '\t', index = False
    )
    # output mi table
    mi_df.to_csv(
        args.output_prefix + '.mi.txt',
        sep = '\t', index = False
    )
    # output removed mismatches
    removed_mismatch_df.to_csv(
        args.output_prefix + '.removed.txt',
        sep = '\t', index = False
    )
    if not variables['mi_calculation_only']:
        # calculate MI empirical p values
        message = 'Score the RNA editing sites.'
        logging.info(message)

        mismatch_df.loc[:, 'mip'] = np.nan

        if mismatch_df['mean_mi'].notna().sum() > 0:
            miecdf = ecdf(
                mismatch_df.loc[
                    mismatch_df['mean_mi'].notna() &
                    (mismatch_df['type'] == 'het_snp'),
                    'mean_mi'
                ]
            )
            mismatch_df.loc[:, 'mip'] = mismatch_df.apply(
                lambda a: miecdf(a['mean_mi'])
                if not np.isnan(a['mean_mi']) else np.nan,
                axis = 1
            )
        else:
            pass

        # run GLM model
        glmresult, scorepf = run_glm(
            mismatch_df,
            mip_threshold = variables['mip_threshold'],
            model = variables['model']
        )

        del glmresult['label']
        del glmresult['id']
        # output mismatch table
        glmresult.to_csv(
            args.output_prefix + '.mismatch.txt',
            sep = '\t', index = False
        )
        scorepf.to_csv(
            args.output_prefix + '.score_performance.txt',
            sep = '\t', index = False
        )
    else:
        pass
    message = 'All done!'
    logging.info(message)


if __name__ == '__main__':
    main()
########################################
