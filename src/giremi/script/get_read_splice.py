#! /usr/bin/env python3
import re
import pysam
import argparse
import multiprocessing as mp
from functools import partial
from collections import defaultdict
from giremi.cs import CS
from giremi.footprint import get_footprints


def get_read_splice_from_bam_in_region(chromosome, start_pos, end_pos,
                                       sam, genome,
                                       mode = 'cs'):
    '''
    Obtain splicing site coordinates from the cs tags in the bam file.
    '''

    splice_list = list()
    for read in sam.fetch(chromosome, start_pos, end_pos):
        read_strand = '-' if read.is_reverse else '+'
        # 0 based
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
            read_seq = read.query_sequence    # reverse complemented for the reversed reads
            ref_seq = genome.fetch(
                chromosome, read.reference_start, read.reference_end
            )
            read_CS = CS.from_cigar_string(
                cigar_string, md_tag_string, read_seq, ref_seq,
                chromosome, read.reference_start, read_strand
            )
        read_introns = read_CS.get_introns(coordinate = 'contig')
        read_introns.sort(key = lambda a: a[0])
        for intron in read_introns:
            splice_list.append(
                [read.query_name,
                 chromosome, intron[0], 'intron_left']
            )
            splice_list.append(
                [read.query_name,
                 chromosome, intron[1], 'intron_right']
            )
    return splice_list


def footprint_bulk_calculation(footprints, variables):
    sam = pysam.AlignmentFile(variables['bam_file'], 'rb')
    if variables['mode'] != 'cs':
        genome = pysam.FastaFile(variables['genome_file'])
    else:
        genome = None

    splice_list = []
    for chromosome, start_pos, end_pos, rc in footprints:
        read_splice = get_read_splice_from_bam_in_region(
            chromosome, start_pos, end_pos,
            sam, genome,
            mode = variables['mode']
        )
        splice_list.extend(read_splice)
    return splice_list


def parse_args():
    parser = argparse.ArgumentParser(
        description="Get read and splice sites from bam file"
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
        default=None
    )
    parser.add_argument(
        "--min_total_depth",
        help = "Min mismatch total read count to be considered (default: 6)",
        type = int,
        default = 2
    )
    parser.add_argument(
        "--mode",
        help = "Mode to find mismatches in BAM: with cs tags (cs), or CIGAR and MD tags (cigar). (defalt: cs)",
        type = str,
        default = 'cs',
        choices = ['cs', 'cigar']
    )
    args = parser.parse_args()
    return args


def main():
    args = parse_args()

    # cigar mode requires the genome fasta.
    if args.mode == 'cigar' and len(args.genome_fasta) < 1:
        raise ValueError('If the mode is "cigar", the genome_fasta must be provided.')
    else:
        pass

    variables = {
        'bam_file': args.bam_file,
        'genome_file': args.genome_fasta,
        'outprefix': args.output_prefix,
        'min_total_depth' : args.min_total_depth,
        'mode': args.mode
    }

    footprints = get_footprints(
        variables['bam_file'], args.chromosomes, variables['min_total_depth']
    )
    n = int(len(footprints) / args.thread)
    footprint_chunk = []
    for i in range(0, len(footprints), n):
        footprint_chunk.append(footprints[i:(i + n)])

    # calculate introns
    with mp.Pool(args.thread) as p:
        footprint_result = p.map(
            partial(footprint_bulk_calculation,
                    variables = variables),
            footprint_chunk
        )

    splicing_list = []
    for spl in footprint_result:
        splicing_list.extend(spl)

    with open('{}.read_splice.txt'.format(variables['outprefix']), 'w') as f:
        f.write(
            '\t'.join([
                'read_name', 'chromosome',
                'pos', 'type'
            ]) + '\n'
        )
        lines = []
        for record in splicing_list:
            lines.append(
                '\t'.join([str(a) for a in record]) +
                '\n'
            )
            if len(lines) > 10000:
                f.writelines(lines)
                lines = []
        f.writelines(lines)


if __name__ == '__main__':
    main()
####################
