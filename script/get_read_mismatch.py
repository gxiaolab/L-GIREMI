#! /usr/bin/env python3
import re
import sys
import pysam
import argparse
import pandas as pd
import numpy as np
from collections import defaultdict

def split_cs_string(cs_string):
    return list(
        zip(
            re.sub('[0-9a-z]', ' ', cs_string).split(),
            re.sub('[:*\-+~]', ' ', cs_string).split()
        )
    )


cslenfuncs = {
    ':': int,
    '*': lambda x: 1,
    '+': lambda x: 0,
    '-': len,
    '~': lambda x: int(re.sub('[a-z]', '', x))
}


def cs_to_df(cs_string, pos):
    cs = split_cs_string(cs_string)
    cslist = list()
    for a, b in cs:
        low = pos
        pos += cslenfuncs[a](b)
        high = pos
        cslist.append([low, high, a, b])
    csdf = pd.DataFrame(
        np.row_stack(cslist),
        columns=['low', 'high', 'ope', 'val']
    )
    csdf.loc[:,'low'] = csdf['low'].astype(int)
    csdf.loc[:,'high'] = csdf['high'].astype(int)
    return csdf


def get_read_mismatch_from_bam(bam_file, chrom,
                               mapq_threshold = 20):
    """
    Obtain mismatch coordinates from the cs tags in the bam file
    """
    mm_dict = defaultdict(dict)
    sam = pysam.AlignmentFile(bam_file, 'rb')
    for read in sam.fetch(chrom):
        if read.mapq < mapq_threshold:
                continue
        elif read.is_secondary: ### skip secondary reads
                continue
        # 0 based
        pos = read.reference_start
        cs = cs_to_df(read.get_tag('cs'), pos)
        cs_splice = cs.loc[cs['ope']=='~']
        cs_mismatch = cs.loc[cs['ope']=='*']
        for ri, row in cs_mismatch.iterrows():
            pos = int(row['low'])
            REF, ALT = row['val'].upper()
            yield [read.query_name, chrom, pos, REF, ALT]


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
    args = parser.parse_args()
    sites = pd.read_csv(
        args.site_file,
        header=None, sep='\t'
    )
    sites.columns = ['chromosome', 'pos']
    sites.loc[:,'sites'] = sites['chromosome'] + ' ' + sites['pos'].astype(str)

    outfile = open('.'.join([args.output_prefix, 'read_mismatch']), 'w')
    outfile.write(
        '\t'.join([
            'read_name', 'chromosome',
            'pos', 'ref', 'alt'
        ]) + '\n'
    )
    for chrom in args.chromosomes:
        lines = []
        aim_pos = sites.loc[sites['chromosome'] == chrom, 'pos'].values
        read_mismatches = get_read_mismatch_from_bam(
            args.bam_file, chrom,
            mapq_threshold = args.mapq_threshold
        )
        for r_m in read_mismatches:
            if r_m[2] in aim_pos:
                lines.append(
                    '\t'.join([str(a) for a in r_m]) +
                    '\n'
                )
            if len(lines) > 10000:
                outfile.writelines(lines)
                lines = []
        outfile.writelines(lines)

    outfile.close()
