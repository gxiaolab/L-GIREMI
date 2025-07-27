from sklearn.metrics import mutual_info_score
from collections import defaultdict
from itertools import combinations


def mismatch_pair_mutual_info(mismatches,
                              min_common_reads=5):
    # mismatches dict should has: ref, type, depth, nt, up, down
    milist = []
    positions = list(mismatches.keys())
    positions.sort()
    for p1, p2 in combinations(positions, 2):
        type1 = mismatches[p1]['type']
        type2 = mismatches[p2]['type']
        p1_read_nt = dict([[rn, nt] for nt, rnlist in mismatches[p1]['nt'].items() for rn in rnlist])
        p2_read_nt = dict([[rn, nt] for nt, rnlist in mismatches[p2]['nt'].items() for rn in rnlist])
        common_read = [rn for rn in p1_read_nt.keys() if rn in p2_read_nt.keys()]
        common_read.sort()
        if len(common_read) < min_common_reads:
            continue
        else:
            p1_nts = [p1_read_nt[rn] for rn in common_read]
            p2_nts = [p2_read_nt[rn] for rn in common_read]
            # convert nts into numbers
            p1_allele_depth = [[nt, mismatches[p1]['depth'][nt]] for nt in mismatches[p1]['depth']]
            p2_allele_depth = [[nt, mismatches[p2]['depth'][nt]] for nt in mismatches[p2]['depth']]
            p1_allele_depth.sort(key = lambda a: a[1], reverse = True)
            p2_allele_depth.sort(key = lambda a: a[1], reverse = True)
            p1_major_allele = p1_allele_depth[0][0]
            p1_minor_allele = p1_allele_depth[1][0]
            p2_major_allele = p2_allele_depth[0][0]
            p2_minor_allele = p2_allele_depth[1][0]
            p1_nt_dict = defaultdict(int)
            p1_nt_dict[p1_minor_allele] = 1
            p1_nt_dict[p1_major_allele] = 2
            p2_nt_dict = defaultdict(int)
            p2_nt_dict[p2_minor_allele] = 1
            p2_nt_dict[p2_major_allele] = 2
            p1_nt_numbers = [p1_nt_dict[nt] for nt in p1_nts]
            p2_nt_numbers = [p2_nt_dict[nt] for nt in p2_nts]
            p1p2mi = mutual_info_score(p1_nt_numbers, p2_nt_numbers)
            milist.append(
                [p1, type1, p2, type2, p1p2mi]
            )
    return milist


def mean_mismatch_pair_mutual_info(mismatch_pair_mi):
    mismatch_mi_lists = defaultdict(list)
    for p1, type1, p2, type2, p1p2mi in mismatch_pair_mi:
        mismatch_mi_lists[p1].append([p2, type2, p1p2mi])
        mismatch_mi_lists[p2].append([p1, type1, p1p2mi])

    mean_mi = defaultdict(int)
    for pos in mismatch_mi_lists:
        mean_mi[pos] = sum(
            [a[2] for a in mismatch_mi_lists[pos]]
        ) / len(mismatch_mi_lists[pos])

    return list([a, b] for a, b in mean_mi.items())

########################################
