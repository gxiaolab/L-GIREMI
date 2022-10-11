from collections import defaultdict


def read_simple_repeat_intervals(repeat_file):
    '''
    Read files with simple repeats
    '''
    f_sr_region = open(repeat_file)
    simple_repeat_intervals = defaultdict(list)
    for line in f_sr_region:
        if line.startswith("#"):
            continue
        line = line.split('\t')
        chrom, start, end = line[0:3]
        start = int(start)
        end = int(end)
        simple_repeat_intervals[chrom].append([start, end])
    f_sr_region.close()
    for chrom in simple_repeat_intervals:
        simple_repeat_intervals[chrom].sort(key = lambda a: a[0])
    return simple_repeat_intervals


def read_snp_positions_in_regions(vcf, chromosome, start_pos, end_pos):
    snp_pos = []
    # location 0 based
    for a in vcf.fetch(chromosome, start_pos, end_pos):
        snp_pos.append(a.start)
    snp_pos.sort()
    return snp_pos


def read_snp_positions_in_chrom(vcf, chromosome):
    snp_pos = []
    # location 0 based
    for a in vcf.fetch(chromosome):
        snp_pos.append(a.start)
    snp_pos.sort()
    return snp_pos

########################################
