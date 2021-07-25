#! /usr/bin/env python3
import re
import sys
import pysam
import argparse
import pandas as pd
import numpy as np
from collections import defaultdict


def get_read_site_from_bam(site, bam_file,
                           mapq_threshold = 20):
    """
    Obtain mismatch coordinates from the cs tags in the bam file
    """
    sam = pysam.AlignmentFile(bam_file, 'rb')
    for ri, row in site.iterrows():
        for pileupcolumn in sam.pileup(
                row['chromosome'],
                row['pos'], row['pos'] + 1,
                truncate=True
        ):
            for read in pileupcolumn.pileups:
                if read.alignment.mapq < mapq_threshold:
                    continue
                elif read.alignment.is_secondary:
                    ### skip secondary reads
                    continue
                if not read.is_del and not read.is_refskip:
                    read_name = read.alignment.query_name
                    read_seq = read.alignment.query_sequence[
                        read.query_position
                    ]
                    if len(read_seq) > 0:
                        yield [
                            read_name,
                            row['chromosome'], row['pos'],
                            read_seq
                        ]


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Get read and mismatch site information from bam file"
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
    args = parser.parse_args()
    sites = pd.read_csv(
        args.site_file,
        header=None, sep='\t'
    )
    sites.columns = ['chromosome', 'pos']

    outfile = open('.'.join([args.output_prefix, 'read_site']), 'w')
    outfile.write(
        '\t'.join([
            'read_name', 'chromosome', 'pos', 'seq'
        ]) + '\n'
    )
    read_sites = get_read_site_from_bam(
        sites, args.bam_file,
        mapq_threshold = args.mapq_threshold
    )
    lines = []
    for r_s in read_sites:
        lines.append('\t'.join([str(a) for a in r_s]) + '\n')
        if len(lines) > 10000:
            outfile.writelines(lines)
            lines = []

    outfile.writelines(lines)
    outfile.close()
