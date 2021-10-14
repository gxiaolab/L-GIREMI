#! /usr/bin/env python3
import argparse
import pandas as pd
from functools import reduce
from collections import defaultdict
from collections import Counter
from sklearn.metrics import mutual_info_score


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description = 'calculate the mutual information of mismatch site and splice site pairs'
    )
    parser.add_argument(
        "-m", "--read_site",
        help = "read-site file: [read_name, chromosome, pos, seq]",
        type=str
    )
    parser.add_argument(
        "-s", "--read_splice",
        help = "corrected read-splice file: [read_name, chromosome, pos, type, corrected_pos, annotation]",
        type=str
    )
    parser.add_argument(
        "-o", "--output_prefix",
        help = "prefix of output file",
        type=str,
        default = 'out'
    )
    args = parser.parse_args()

    site_file = args.read_site
    splice_file = args.read_splice
    pair_out_file = args.output_prefix + '.site_splice_pair'

    # count reads covers site-splice pairs
    rsite = pd.read_table(site_file, header = 0, sep = '\t')

    rsplice = pd.read_table(
        splice_file, header = 0, sep = '\t',
        chunksize = 10000
    )

    site_splice_pair = Counter()

    # calculation by chunk
    for ckdata in rsplice:
        mdata = pd.merge(
            rsite,
            ckdata[['read_name', 'chromosome', 'corrected_pos']],
            how = 'inner', on = ['read_name', 'chromosome']
        )
        mdata.columns = ['read_name', 'chromosome', 'site_pos', 'seq', 'splice_pos']

        counts = mdata.groupby(
            ['chromosome', 'site_pos', 'seq', 'splice_pos']
        )['read_name'].count().reset_index()

        for ri, row in counts.iterrows():
            pid = '\t'.join(
                [
                    str(a)
                    for a in list(
                            row[['chromosome', 'site_pos', 'seq',
                                 'splice_pos']])
                ]
            )
            site_splice_pair[pid] += row['read_name']

    tmppairs = pd.DataFrame(
        {'key': key, 'count': val}
        for key, val in site_splice_pair.items()
    )

    pairs = tmppairs['key'].str.split('\t', expand = True)

    pairs.columns = [
        'chromosome', 'site_pos', 'seq', 'splice_pos'
    ]

    pairs.loc[:, 'count'] = tmppairs['count']

    # site dict
    sites = defaultdict(lambda: defaultdict(list))

    for ri, row in rsite.iterrows():
        poslabel = '{}:{}'.format(row['chromosome'], row['pos'])
        sites[poslabel][row['seq']].append(row['read_name'])

    # splice dict
    splices = defaultdict(list)

    with open(splice_file, 'r') as f:
        for line in f:
            cols = line.strip().split('\t')
            poslabel = ':'.join([cols[1], cols[4]])
            splices[poslabel].append(cols[0])

    pairs.loc[:, 'mi'] = -1

    for ri, row in pairs.iterrows():
        # 'chromosome'
        # 'site_pos'
        # 'seq'
        # 'splice_pos'
        site_label = ':'.join(
            [row['chromosome'], str(row['site_pos'])]
        )
        splice_label = ':'.join(
            [row['chromosome'], str(row['splice_pos'])]
        )
        read_site_seq = sites[site_label][row['seq']]
        read_splice = splices[splice_label]
        read_site_all = sorted(reduce(
            lambda a, b: a + b,
            [sites[site_label][c] for c in sites[site_label]]
        ))
        with_seq = [int(a in read_site_seq) for a in read_site_all]
        with_splice = [int(a in read_splice) for a in read_site_all]
        pairs.loc[ri, 'mi'] = mutual_info_score(with_seq, with_splice)

    # output
    pairs.to_csv(
        pair_out_file, sep = '\t', index = False
    )

####################
