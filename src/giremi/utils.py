import numpy as np


def merge_intervals(intervals):
    sorted_intervals = intervals.copy()
    sorted_intervals.sort(key = lambda a: a[0])
    merged = list()
    interval_counts = list()
    for interval in sorted_intervals:
        if not merged or merged[-1][1] < interval[0]:
            merged.append(interval)
            interval_counts.append(1)
        else:
            merged[-1][1] = max(merged[-1][1], interval[1])
            interval_counts[-1] += 1
    return merged, interval_counts


def positions_in_intervals(positions, intervals):
    # check whether the positions in the intervals
    intervals.sort(key = lambda a: a[0])
    intervals_start = [a for a, b in intervals]
    intervals_end = [b for a, b in intervals]
    # find the start and end insertion position for positions
    start_idx = np.searchsorted(intervals_start, positions, side = 'right')
    end_idx = np.searchsorted(intervals_end, positions, side = 'right')
    # only positions with start - end = 1 in some of the intervals
    in_intervals = [a - b == 1 for a, b in zip(start_idx, end_idx)]
    return in_intervals, [b if a else -1 for a, b in zip(in_intervals, end_idx)]

########################################
