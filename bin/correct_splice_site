#! /usr/bin/env python3
import pysam
import argparse
import numpy as np
from collections import defaultdict


def get_gtf_splice_pos(gtf_file, chromosomes):
    gtf = pysam.TabixFile(gtf_file)

    pos = dict()

    poslist = defaultdict(list)

    for chrom in chromosomes:
        for gtf_entry in gtf.fetch(
                chrom,
                parser=pysam.asGTF()):
            if gtf_entry.feature == 'exon':
                poslist[chrom].append([gtf_entry.start, 'exon_start'])
                poslist[chrom].append([gtf_entry.end, 'exon_end'])
        pos[chrom] = np.sort(np.unique(np.array(
            [a for a, b in poslist[chrom] + poslist[chrom]]
        )))

    return pos


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Correct splice sites by gtf file"
    )
    parser.add_argument(
        "-s", "--splice_file",
        help="input splice file [read_name, chromosome, pos, type] with column names",
        type=str,
        default=None,
        required=True
    )
    parser.add_argument(
        "--annotation_gtf",
        help = "gtf (gz and tabix indexed) file of genome annotation (gencode)",
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
        "--window",
        help = "window to be considered as the same splice site (default: 10)",
        type=int,
        default=10
    )
    args = parser.parse_args()

    pos = get_gtf_splice_pos(args.annotation_gtf,
                             args.chromosomes)

    outfile = open('.'.join([args.output_prefix, 'corrected_read_splice']), 'w')
    outfile.write(
        '\t'.join([
            'read_name', 'chromosome',
            'pos', 'type', 'corrected_pos', 'annotation'
        ]) + '\n'
    )
    outlines = []
    with open(args.splice_file, 'r') as f:
        for theline in f:
            if theline.find('read_name') >= 0:
                continue
            cols = theline.strip().split('\t')
            thechrom = cols[1]
            thepos = int(cols[2])
            thepostype = cols[3]
            idx = abs(pos[thechrom] - thepos) <= args.window
            idx_true = np.ravel(np.argwhere(idx))
            corrected_pos = thepos
            ingtf = 'FALSE'
            if len(idx_true) > 0:
                corrected_pos = pos[thechrom][idx_true[0]]
                ingtf = 'TRUE'
            outlines.append(
                '\t'.join(cols + [str(corrected_pos), ingtf]) +
                '\n'
            )
            if len(outlines) > 10000:
                outfile.writelines(outlines)
                outlines = []
        outfile.writelines(outlines)

    outfile.close()

####################
