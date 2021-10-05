#! /usr/bin/env python3
import re
import pysam
import argparse
import pandas as pd
import numpy as np
import multiprocessing as mp
from functools import partial


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


def chrom_get_read_intron_from_bam(chrom, variables):
    """
    Obtain intron coordinates from the cs tags in the bam file
    """
    intron_list = list()
    sam = pysam.AlignmentFile(variables['bam_file'], 'rb')
    for read in sam.fetch(chrom):
        if read.mapq < variables['mapq_threshold']:
            continue
        elif read.is_secondary:
            # skip secondary reads
            continue
        # 0 based
        pos = read.reference_start
        cs = cs_to_df(read.get_tag('cs'), pos)
        cs_splice = cs.loc[cs['ope'] == '~']
        for ri, row in cs_splice.iterrows():
            low = int(row['low'])
            high = int(row['high'])
            a = row['val'][0:2]
            b = row['val'][-2:]
            length = int(row['val'][2:-2])
            intron_list.append(
                [read.query_name,
                 chrom, low, high, length,
                 a, b]
            )
    return intron_list


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Get read and intron information from bam file"
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
        'outprefix': args.output_prefix,
        'mapq_threshold': args.mapq_threshold
    }

    with mp.Pool(args.thread) as p:
        results = p.map(
            partial(
                chrom_get_read_intron_from_bam,
                variables=variables
            ),
            args.chromosomes
        )

    outfile = open('.'.join([variables['outprefix'], 'read_intron']), 'w')
    outfile.write(
        '\t'.join([
            'read_name', 'chromosome',
            'start', 'end', 'length', 'a', 'b'
        ]) + '\n'
    )
    for chrom_result in results:
        lines = []
        for r_m in chrom_result:
            lines.append(
                '\t'.join([str(a) for a in r_m]) +
                '\n'
            )
            if len(lines) > 10000:
                outfile.writelines(lines)
                lines = []
        outfile.writelines(lines)

    outfile.close()

####################
