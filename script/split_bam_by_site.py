#! /usr/bin/env python3
import pysam
import argparse
from collections import defaultdict


def get_seq_dict(sam, chrom, pos):

    seqdict = defaultdict(list)

    for pileupcolumn in sam.pileup(contig = chrom,
                                   start = pos,
                                   end = pos + 1):
        if pileupcolumn.pos == pos:
            read_names = pileupcolumn.get_query_names()
            read_nts = [
                a.upper() for a in pileupcolumn.get_query_sequences()
            ]
            nt_set = list(set(read_nts))
            idx = list(range(len(read_names)))
            for nt in nt_set:
                nt_idx = [
                    i for i in idx
                    if read_nts[i] == nt
                ]
                seqdict[nt] = [
                    pileupcolumn.pileups[i].alignment
                    for i in nt_idx
                ]
    return seqdict

####################


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Save reads that cover one position into separated SAM files by the genomic location"
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
        "-c", "--chromosome",
        type = str
    )
    parser.add_argument(
        "-p", "--pos",
        type = int
    )
    args = parser.parse_args()

    sam = pysam.AlignmentFile(args.bam_file, 'rb')
    seqdict = get_seq_dict(
        sam, args.chromosome, args.pos
    )

    for nt in seqdict.keys():
        outfile = pysam.AlignmentFile(
            args.output_prefix + "_" + nt + '.sam',
            "w", template=sam
        )
        for s in seqdict[nt]:
            outfile.write(s)
        outfile.close()
    sam.close()

####################
