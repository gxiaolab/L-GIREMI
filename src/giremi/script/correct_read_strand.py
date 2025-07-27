#! /usr/bin/env python3
import argparse
import pysam
import logging
import pandas as pd
import multiprocessing as mp
from functools import partial
from giremi.footprint import get_footprints
from giremi.strand import correct_read_strand_in_region

def footprint_bulk_calculation(footprints, variables):
    sam = pysam.AlignmentFile(variables['bam_file'], 'rb')
    genome = pysam.FastaFile(variables['genome_file'])
    gtf = pysam.TabixFile(variables['gtf_file'])

    strand_list = []
    for footprint in footprints:
        chromosome, start_pos, end_pos, rc = footprint
        read_strand_list = correct_read_strand_in_region(
            chromosome, start_pos, end_pos,
            sam, gtf, genome,
            variables['gene_padding'], variables['exon_padding'],
            keep_non_spliced_read = variables['keep_non_spliced_read'],
            mode = variables['mode']
        )
        strand_list.extend(read_strand_list)
    strand_df = pd.DataFrame.from_records(
        strand_list,
        columns = ['read_name', 'original_read_strand', 'corrected_read_strand']
    )
    sam.close()
    genome.close()
    gtf.close()
    return strand_df


def parse_args():
    parser = argparse.ArgumentParser(
        description="Correct read strand."
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
        type = str,
        default = 'out'
    )
    parser.add_argument(
        "-t", "--thread",
        help = "cores to be used",
        type = int,
        default = 1
    )
    parser.add_argument(
        "--annotation_gtf",
        help = "gtf (gz and tabix indexed) file of genome annotation (gencode)",
        type = str,
        default = None,
        required = True
    )
    parser.add_argument(
        "--genome_fasta",
        help = "path of genome fasta file",
        type = str,
        default = None,
        required = True
    )
    parser.add_argument(
        "--keep_non_spliced_read",
        help = "Keep non spliced reads (default: without this parameter, reads without splicing would be dropped.)",
        action = 'store_true'
    )
    parser.add_argument(
        "--mode",
        help = "Mode to find mismatches in BAM: with cs tags (cs), or CIGAR and MD tags (cigar). (defalt: cs)",
        type = str,
        default = 'cs'
    )
    parser.add_argument(
        "--padding_exon",
        help = "expand the range when searching exon gtf (default: 10)",
        type = int,
        default = 10
    )
    parser.add_argument(
        "--padding_gene",
        help = "expand the range when searching gene gtf (default: 500)",
        type = int,
        default = 500
    )
    args = parser.parse_args()
    return args


def main():
    args = parse_args()
    variables = {
        'bam_file' : args.bam_file,
        'genome_file' : args.genome_fasta,
        'gtf_file' : args.annotation_gtf,
        'exon_padding' : args.padding_exon,
        'gene_padding' : args.padding_gene,
        'keep_non_spliced_read' : args.keep_non_spliced_read,
        'mode': args.mode
    }

    # logging

    logging.basicConfig(
        encoding='utf-8',
        format='%(asctime)s %(levelname)s %(message)s',
        level=logging.INFO
    )

    # footprint depth use 1
    message = 'Get regions with reads.'
    logging.info(message)
    footprints = get_footprints(
        variables['bam_file'], args.chromosomes, 1
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
                    variables = variables),
            footprint_chunk
        )
    strand_df_list = []
    for dfstrand in footprint_result:
        strand_df_list.append(dfstrand)
    strand_df = pd.concat(strand_df_list, axis = 0)

    # output strand table
    strand_df.to_csv(
        args.output_prefix + '.strand.txt',
        sep = '\t', index = False
    )
    message = 'All done!'
    logging.info(message)


if __name__ == '__main__':
    main()
########################################
