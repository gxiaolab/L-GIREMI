#! /usr/bin/env python3
import pysam
import argparse
import pandas as pd
import multiprocessing as mp
from functools import partial


def chrom_get_read_site_from_bam(chrom, variables):
    """
    Obtain mismatch coordinates from the cs tags in the bam file
    """
    sites = pd.read_csv(
        variables['aim_site_file'],
        header=None, sep='\t'
    ).drop_duplicates()
    sites.columns = ['chromosome', 'pos']
    sites = sites.loc[sites['chromosome'] == chrom]
    rs_list = list()
    sam = pysam.AlignmentFile(variables['bam_file'], 'rb')
    for ri, row in sites.iterrows():
        for pileupcolumn in sam.pileup(
                row['chromosome'],
                row['pos'], row['pos'] + 1,
                truncate=True
        ):
            for read in pileupcolumn.pileups:
                if read.alignment.mapq < variables['mapq_threshold']:
                    continue
                elif read.alignment.is_secondary:
                    # skip secondary reads
                    continue
                if not read.is_del and not read.is_refskip:
                    read_name = read.alignment.query_name
                    read_seq = read.alignment.query_sequence[
                        read.query_position
                    ]
                    if len(read_seq) > 0:
                        rs_list.append(
                            [read_name,
                             row['chromosome'],
                             row['pos'],
                             read_seq]
                        )
    return rs_list


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Get read and site information from bam file"
    )
    parser.add_argument(
        "-s", "--site_file",
        help="input site file: [chromosome, pos], without header, separated by tab",
        type=str,
        default=None,
        required=True
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
        "--mapq_threshold",
        help = "Min MAPQ to be considered in bam file (default: 20)",
        type=float,
        default=20
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
        'aim_site_file': args.site_file,
        'outprefix': args.output_prefix,
        'mapq_threshold': args.mapq_threshold
    }

    with mp.Pool(args.thread) as p:
        results = p.map(
            partial(
                chrom_get_read_site_from_bam,
                variables=variables
            ),
            args.chromosomes
        )

    outfile = open('.'.join([variables['outprefix'], 'read_site']), 'w')
    outfile.write(
        '\t'.join([
            'read_name', 'chromosome', 'pos', 'seq'
        ]) + '\n'
    )
    for chrom_result in results:
        lines = []
        for r_s in chrom_result:
            lines.append('\t'.join([str(a) for a in r_s]) + '\n')
            if len(lines) > 10000:
                outfile.writelines(lines)
                lines = []
        outfile.writelines(lines)
    outfile.close()

####################
