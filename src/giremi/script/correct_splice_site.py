#! /usr/bin/env python3
import pysam
import argparse
import multiprocessing as mp
from functools import partial


def get_gtf_splice_pos_in_region(gtf, chromosome, start_pos, end_pos):
    poslist = list()

    for gtf_entry in gtf.fetch(
            chromosome, start_pos, end_pos,
            parser = pysam.asGTF()
    ):
        if gtf_entry.feature == 'exon' or gtf_entry.feature == 'CDS':
            poslist.append([gtf_entry.start, 'exon_start'])
            poslist.append([gtf_entry.end, 'exon_end'])

    return poslist


def bulk_calculation(chunk, variables):
    gtf = pysam.TabixFile(variables['gtf_file'])
    splicing_info_list = []
    for rname, chrom, pos, postype in chunk:
        gtf_pos_list = get_gtf_splice_pos_in_region(
            gtf, chrom,
            pos - variables['window'],
            pos + variables['window']
        )
        if len(gtf_pos_list) == 0:
            splicing_info_list.append(
                [rname, chrom, pos, postype,
                 -1, '']
            )
        else:
            dist = [[abs(gpos - pos), gpos, gtype]
                    for gpos, gtype in gtf_pos_list]
            dist.sort(key = lambda a: a[0])
            if dist[0][0] <= variables['window']:
                splicing_info_list.append(
                    [rname, chrom, pos, postype,
                     dist[0][1], dist[0][2]]
                )
            else:
                splicing_info_list.append(
                    [rname, chrom, pos, postype,
                     -1, '']
                )

    gtf.close()
    return splicing_info_list


def parse_args():
    parser = argparse.ArgumentParser(
        description="Correct splice sites by gtf file"
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
        "-s", "--splice_file",
        help="input splice file [read_name, chromosome, pos, type] with column names",
        type=str,
        default=None,
        required=True
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
        type=str,
        default=None,
        required=True
    )
    parser.add_argument(
        "--window",
        help = "window to be considered as the same splice site (default: 10)",
        type=int,
        default=10
    )
    args = parser.parse_args()
    return args


def main():
    args = parse_args()

    variables = {
        'gtf_file' : args.annotation_gtf,
        'window': args.window
    }

    read_splicing = []
    with open(args.splice_file, 'r') as f:
        for theline in f:
            if theline.find('read_name') >= 0:
                continue
            else:
                rname, chrom, pos, postype = theline.strip().split('\t')
                read_splicing.append([rname, chrom, int(pos), postype])
    n = int(len(read_splicing) / args.thread / 2)
    read_splicing_chunk = []
    for i in range(0, len(read_splicing), n):
        read_splicing_chunk.append(read_splicing[i: (i + n)])

    with mp.Pool(args.thread) as p:
        chunk_result = p.map(
            partial(bulk_calculation,
                    variables = variables),
            read_splicing_chunk
        )

    corrected_splicing_list = list()
    for result in chunk_result:
        corrected_splicing_list.extend(result)

    with open('{}.corrected_read_splice.txt'.format(args.output_prefix), 'w') as f:
        f.write(
            '\t'.join([
                'read_name', 'chromosome',
                'pos', 'type', 'corrected_pos', 'annotation'
            ]) + '\n'
        )
        lines = []
        for record in corrected_splicing_list:
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
