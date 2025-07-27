# Here, we define the regions that overlapped reads covered as footprints.
import pysam
from giremi.utils import merge_intervals


def chr_get_footprints(sam, chromosome, min_read_count = 2):
    '''
    Get the footprints of reads (aka overlapped regions with reads) from one chromosome.

    Parameters:
        sam            - pysam AlignmentFile object.
        chromosome     - string, name of the chromosome to process.
        min_read_count - integer, minimum read counts for footprints to return.

    Return: (footprints, read_counts)
        footprints  - list of footprints intervals (start, end) that passed the min_read_count threshold.
        read_counts - read counts for the intervals that passed the min_read_count threshold.
    '''
    read_intervals = []
    for read in sam.fetch(chromosome):
        # 0 based
        read_start = read.reference_start
        read_end = read.reference_end
        read_intervals.append([read_start, read_end])
    footprints, read_counts = merge_intervals(read_intervals)
    idx = [a >= min_read_count for a in read_counts]
    return [f for f, a in zip(footprints, idx) if a], \
        [r for r, a in zip(read_counts, idx) if a]


def get_footprints(bam_file, chromosomes, min_read_count = 2):
    '''
    Get all the footprints from the chromosomes listed.

    Parameters:
        bam_file       - string, file path for the bam file to process.
        chromosome     - string, name of the chromosome to process.
        min_read_count - integer, minimum read counts for footprints to return.

    Return: list of footprint information (each footprint has a list: [chromosome, start, end, read_counts]).
    '''
    sam = pysam.AlignmentFile(bam_file, 'rb')
    intervals = []
    for chrom in chromosomes:
        footprints, read_counts = chr_get_footprints(sam, chrom, min_read_count)
        intervals.extend(
            [[chrom, f[0], f[1], r] for f, r in zip(footprints, read_counts)]
        )
    sam.close()
    return intervals


########################################
