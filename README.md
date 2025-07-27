# L-GIREMI

[![](https://img.shields.io/badge/version-v0.2.4-blue)](https://pypi.org/project/l-giremi/)

L-GIREMI (Long-read Genome-independent Identification of RNA Editing by Mutual Information) is a method for identification of RNA editing sites from long-read RNA-seq data.

## Requirement

The L-GIREMI software was developed with python3 on Linux, which
demands several python packages.

* python3.5+: with sys, argparse, re, logging, functools, collections,
  multiprocessing.
* [scikit-learn](https://scikit-learn.org): 0.20+
* [scipy](https://www.scipy.org): 1.5+
* [numpy](https://numpy.org): 1.10+
* [pandas](https://pandas.pydata.org): 1.0+
* [pysam](https://pysam.readthedocs.io): 0.16+

And the analysis with L-GIREMI also needs additional software:
* [samtools](http://www.htslib.org)
* [bcftools](http://www.htslib.org)
* [bgzip](http://www.htslib.org)
* [tabix](http://www.htslib.org)
* [minimap2](https://github.com/lh3/minimap2)

### Installation

It's recommended to install the software into a virtual environment.

Create virtual environment:

```{bash}
conda create -n l_giremi
conda activate l_giremi
```

Or:

```{bash}
virtualenv l_giremi
source bin/activate
```


#### From github

Download from github:

```{bash}
git clone https://github.com/gxiaolab/L-GIREMI
cd L-GIREMI
```

Install the package:

```{bash}
pip install .
```

### Reference data

Users should make sure to use reference files with the same assembly
version and the same chromosome naming convention in all the
processes.

#### Genome fasta file

L-GIREMI requires a genome fasta file as input. For the human genome,
the fasta file can be obtained from
[UCSC](https://hgdownload.cse.ucsc.edu/goldenpath/hg38/chromosomes/)
or [NCBI](https://www.ncbi.nlm.nih.gov/genome/guide/human/). The fasta
file should include sequences from all the chromosomes.


#### SNP vcf file

L-GIREMI will use a SNP vcf file to get putative heterozygous SNPs. A
vcf of dbSNP is acceptable. Alternatively, the users can provide a
sample-specific SNP vcf file, if available. The vcf file should be
converted into the bcf format, sorted, and indexed.

The following is an example flow to obtain dbSNP vcf and convert the
dbSNP chromosome names into UCSC format.


```{bash}
wget "https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.40.gz"

wget "https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.40.gz.md5"

md5sum -c GCF_000001405.40.gz.md5

wget "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.28_GRCh38.p13/GCA_000001405.28_GRCh38.p13_assembly_report.txt"

grep -E -v "^#" GCA_000001405.28_GRCh38.p13_assembly_report.txt | cut -f 7,10 | tr "\t" " " > GRCh38-to-hg38.map

bcftools annotate --rename-chrs GRCh38-to-hg38.map GCF_000001405.40.gz -Ob -o dbsnp.40.hg38.bcf

bcftools index dbsnp.40.hg38.bcf
```

#### Gene annotation gtf

A gene annotation gtf is required for L-GIREMI, which can be obtained
from [gencode](https://www.gencodegenes.org/human/). The gtf file
should be sorted, zipped, and indexed.


The following is an example flow to prepare the gtf file:

```{bash}
(zgrep ^"#" gencode.v48.annotation.gtf.gz; \
    zgrep -v ^"#" gencode.v48.annotation.gtf.gz | \
    sort -k1,1 -k4,4n) | \
    bgzip > gencode.v48.annotation.sorted.gtf.gz

tabix -p gff gencode.v48.annotation.sorted.gtf.gz
```

#### Simple repeat region table


A table containing annotated simple repeats is used by L-GIREMI to
pre-filter sites. It can be obtained from [UCSC table
browser](http://genome.ucsc.edu/cgi-bin/hgTables). The table format
should be: chromosome, start, end. No header is needed for the file.

```{bash}
zcat ucsc_hgtables_simplerepeats.tsv.gz |cut -f 2,3,4 |grep -v chrom > repeats.txt
```

## Usage

```{bash}
usage: l-giremi [-h] -b BAM_FILE [-c [CHROMOSOMES ...]] [-o OUTPUT_PREFIX] [-t THREAD] --annotation_gtf ANNOTATION_GTF --genome_fasta GENOME_FASTA [--homopoly_length HOMOPOLY_LENGTH] [--keep_non_spliced_read]
                [--max_het_snp_ratio MAX_HET_SNP_RATIO] [--mi_calculation_only] [--mi_min_common_read MI_MIN_COMMON_READ] [--mi_p_threshold MI_P_THRESHOLD] [--min_allele_depth MIN_ALLELE_DEPTH]
                [--min_allele_ratio MIN_ALLELE_RATIO] [--min_het_snp_ratio MIN_HET_SNP_RATIO] [--min_total_depth MIN_TOTAL_DEPTH] [--min_dist_from_splice MIN_DIST_FROM_SPLICE] [--mode MODE] [--model MODEL]
                [--padding_exon PADDING_EXON] [--padding_gene PADDING_GENE] --repeat_txt REPEAT_TXT [--skip_strand_correction] --snp_bcf SNP_BCF [--version]

L-GIREMI (Long-read RNA-seq Genome-independent Identification of RNA Editing by Mutual Information)

optional arguments:
  -h, --help            show this help message and exit
  -b BAM_FILE, --bam_file BAM_FILE
                        input bam file, with cs tags, sorted and indexed
  -c [CHROMOSOMES ...], --chromosomes [CHROMOSOMES ...]
                        chromosomes to be analyzed
  -o OUTPUT_PREFIX, --output_prefix OUTPUT_PREFIX
                        prefix of output file
  -t THREAD, --thread THREAD
                        cores to be used
  --annotation_gtf ANNOTATION_GTF
                        gtf (gz and tabix indexed) file of genome annotation (gencode)
  --genome_fasta GENOME_FASTA
                        path of genome fasta file
  --homopoly_length HOMOPOLY_LENGTH
                        left and right sequence length to be searched for the homopoly around sites (default: 5)
  --keep_non_spliced_read
                        Keep non spliced reads (default: without this parameter, reads without splicing would be dropped.)
  --max_het_snp_ratio MAX_HET_SNP_RATIO
                        Max ratio to be considered as heterogenous SNPs (default: 0.65)
  --mi_calculation_only
                        Only calculate the mutual information, without the modeling step. (default: calculate the steps after the MI step)
  --mi_min_common_read MI_MIN_COMMON_READ
                        Min common read for site pairs to calculate MI (default: 6)
  --mi_p_threshold MI_P_THRESHOLD
                        MI p value threshold to be used to separate RNA editing sites (default: 0.05)
  --min_allele_depth MIN_ALLELE_DEPTH
                        Min mismatch allele read count to be considered (default: 3)
  --min_allele_ratio MIN_ALLELE_RATIO
                        Min mismatch ratio to be considered (default: 0.05)
  --min_het_snp_ratio MIN_HET_SNP_RATIO
                        Min ratio to be considered as heterogenous SNPs (default: 0.35)
  --min_total_depth MIN_TOTAL_DEPTH
                        Min mismatch total read count to be considered (default: 6)
  --min_dist_from_splice MIN_DIST_FROM_SPLICE
                        Drop sites within the distance from splice junctions (default: 4)
  --mode MODE           Mode to find mismatches in BAM: with cs tags (cs), or CIGAR and MD tags (cigar). (defalt: cs)
  --model MODEL         GLM model for scoring: lm, logistic. (defalt: lm)
  --padding_exon PADDING_EXON
                        expand the range when searching exon gtf (default: 10)
  --padding_gene PADDING_GENE
                        expand the range when searching gene gtf (default: 500)
  --repeat_txt REPEAT_TXT
                        path of txt file of simple repeats [chromosom, start, end] (0 based)
  --skip_strand_correction
                        skip the strand correction step to save runtime.
  --snp_bcf SNP_BCF     path of dbSNP bcf file
  --version             show program's version number and exit
```

## Analysis process

Once the following files are available (see above):
* reference FASTA file: for example hg38 fasta.
* SNP VCF file: dbSNP VCF file (the chromosome names should be agreed
  with the reference FASTA file), or known SNP VCF file for the
  sample.
* genome annotation GTF file: with exon annotation and gene
  annotation, for example GENCODE gtf file.

First, align the long-read RNA-seq fastq file using
[minimap2](https://github.com/lh3/minimap2) with cs tag, which is
required by `l-giremi`.

```{bash}
minimap2 -t 8 -ax splice -uf \
         --secondary=no -N 5 --cs \
         -o SAM_FILE GENOME_FILE FASTQ_FILE
```

Then, remove the unmapped and non-uniquely mapped reads from the SAM
file, and sort it into an indexed BAM file.

```{bash}
samtools view -O BAM -F 2052 -h $SAM_FILE | \
    samtools sort -O BAM -@ 7 -o $BAM_FILE -

samtools index $BAM_FILE
```

Next, run `l-giremi`.

```{bash}
l-giremi \
    -t 8 \
    --bam_file $BAM_FILE \
    --output_prefix $OUTPREFIX \
    --genome_fasta $GENOME_FILE \
    --snp_bcf $SNP_FILE \
    --repeat_txt $REPEAT_FILE \
    --annotation_gtf $GTF_FILE \
    --chromosomes chr1 chr2 chr3 chr4 \
    chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 \
    chr13 chr14 chr15 chr16 chr17 chr18 chr19 \
    chr20 chr21 chr22 chrX chrY
```

### Only MI calculation

If only the MI calculation is needed, for the purpose of saving time,
you can set the parameter to ignore the steps after the MI step.

```{bash}
l-giremi \
    -t 8 \
    --bam_file $BAM_FILE \
    --output_prefix $OUTPREFIX \
    --genome_fasta $GENOME_FILE \
    --snp_bcf $SNP_FILE \
    --repeat_txt $REPEAT_FILE \
    --annotation_gtf $GTF_FILE \
    --chromosomes chr1 chr2 chr3 chr4 \
    chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 \
    chr13 chr14 chr15 chr16 chr17 chr18 chr19 \
    chr20 chr21 chr22 chrX chrY \
    --mi_calculation_only
```

And, only the mutual information result would be output.

## Output

`l-giremi` generates several output files:

* corrected read strand files: the files are stored as `$OUTPREFIX.strand.txt`.
  columns:
  1. read_name: the read name.
  2. corrected\_read\_strand: corrected read strand.
* mutual information files: stored as `$OUTPREFIX.mi.txt`.
  columns:
  1. chromosome: chromosome name.
  2. strand: strand of the mismatch sites.
  3. site1\_pos: position of site 1 on the chromosome, 0-based.
  4. site1\_type: type of site 1: snp, het\_snp, or mismatch. snp and
     het\_snp for sites overlapping the dbSNP annotation, and het\_snp
     site ratio satisfying the parameters.
  5. site2\_pos: position of site 2 on the chromosome, 0-based.
  6. site2\_type: type of site 2: snp, het\_snp, or mismatch.
  7. mi: mutual information value.
* mismatch score file: stored as `$OUTPREFIX.mismatch.txt`.
  columns:
  1. type: het_snp, for mismatch sites that overlapping with SNP
     annotations and with mismatch ratio satisfying the
     parameters. mimatch for other types.
  2. chromosome: chromosome name.
  3. pos: position on the chromosome, 0-based.
  4. strand: strand of the mismatch sites.
  5. change\_type: the mismatch type, [ref]>[alt].
  6. ratio: the mismatch ratio.
  7. allelic\_ratio\_diff: difference between the major mismatch ratio from
     regional allelic ratio (ratio - allelic\_ratio).
  8. depth: the total read count for the site.
  9. A:C:T:G: read counts for A, C, T, G, not strand specific.
  10. up\_seq: a 5' nucleotide ahead of the mismatch sites.
  11. down_seq: a 3' nucleotide after the mismatch sites.
  12. mean\_mi: mean mutual information with the het\_snps, empty if not applicable.
  13. mip: emperical p value of the MI, empty if not applicable.
  14. score: RNA editing score by the model.
* score performance file: stored as `$OUTPREFIX.score_performance`.
  columns:
  1. score: score from mismatch score file
  2. tp: true positive count, the number of sites with score larger
     than or equal to the score in the line used as positive training
     data
  3. fp: false positive count, the number of sites with score larger
     than or equal to the score in the line used as negative training
     data
  4. fn: false negative count, the number of sites with score smaller
     than the score in the line used as positive training data
  5. tn: true negative count, the number of sites with score smaller
     than the score in the line used as negative training data
  6. precision: tp / (tp + fp)
  7. recall: tp / (tp + fn)
  8. f1: 2 * (precision * recall) / (precision + recall)
  9. sensitivity: tp / (tp + fn)
  10. specificity: tn / (tn + fp)
  11. tpr: true positive rate, tp / (tp + fn)
  12. fpr: false positive rate, fp / (fp + tn)
  13. max\_f1: bool, indicating the max f1

## Useful scripts

There are also some useful scripts provided with the main l-giremi
code.

### `correct_read_strand`

`correct_read_strand` is used for the correction of strand of long
reads.  It will use the information of splicing pattern of reads, gtf
annotations to correct the strand. The output is a txt file with 3
columns: read\_name, original\_read\_strand, corrected\_read\_strand.

```{bash}
usage: correct_read_strand [-h] -b BAM_FILE [-c [CHROMOSOMES ...]] [-o OUTPUT_PREFIX] [-t THREAD] --annotation_gtf ANNOTATION_GTF --genome_fasta GENOME_FASTA [--keep_non_spliced_read] [--mode MODE] [--padding_exon PADDING_EXON]
                           [--padding_gene PADDING_GENE]

Correct read strand.

optional arguments:
  -h, --help            show this help message and exit
  -b BAM_FILE, --bam_file BAM_FILE
                        input bam file, with cs tags, sorted and indexed
  -c [CHROMOSOMES ...], --chromosomes [CHROMOSOMES ...]
                        chromosomes to be analyzed
  -o OUTPUT_PREFIX, --output_prefix OUTPUT_PREFIX
                        prefix of output file
  -t THREAD, --thread THREAD
                        cores to be used
  --annotation_gtf ANNOTATION_GTF
                        gtf (gz and tabix indexed) file of genome annotation (gencode)
  --genome_fasta GENOME_FASTA
                        path of genome fasta file
  --keep_non_spliced_read
                        Keep non spliced reads (default: without this parameter, reads without splicing would be dropped.)
  --mode MODE           Mode to find mismatches in BAM: with cs tags (cs), or CIGAR and MD tags (cigar). (defalt: cs)
  --padding_exon PADDING_EXON
                        expand the range when searching exon gtf (default: 10)
  --padding_gene PADDING_GENE
                        expand the range when searching gene gtf (default: 500)
```

### `correct_splice_site`

Since the raw splice sites detected in the long-read RNA-seq data may
not be accurate, a correction procedure is needed [(Wyman and
Mortazavi, 2019)][2]. `correct_splice_site` adopted similar correction
strategies with TranscriptClean. The correction needs a GTF annotation
file.


```{bash}
usage: correct_splice_site [-h] -s SPLICE_FILE --annotation_gtf ANNOTATION_GTF
                           [-c [CHROMOSOMES ...]] [-o OUTPUT_PREFIX] [--window WINDOW]

Correct splice sites by gtf file

optional arguments:
  -h, --help            show this help message and exit
  -s SPLICE_FILE, --splice_file SPLICE_FILE
                        input splice file [read_name, chromosome, pos, type] with column names
  --annotation_gtf ANNOTATION_GTF
                        gtf (gz and tabix indexed) file of genome annotation (gencode)
  -c [CHROMOSOMES ...], --chromosomes [CHROMOSOMES ...]
                        chromosomes to be analyzed
  -o OUTPUT_PREFIX, --output_prefix OUTPUT_PREFIX
                        prefix of output file
  --window WINDOW       window to be considered as the same splice site (default: 10)
```


### `get_aei`

`get_aei` is used for the calculation of Alu Editing Index (AEI)
from long-read RNA-seq data [(Roth et al., 2019)][1]. Calculation of
AEI neads the read strand files (output of `l-giremi`), the genome
fasta file, a txt file with Alu locations, and the SNP BCF reference.

```{bash}
usage: get_aei [-h] -b BAM_FILE --strand_file STRAND_FILE [-c [CHROMOSOMES ...]] [-o OUTPUT_PREFIX]
                  --genome_fasta GENOME_FASTA --snp_bcf SNP_BCF --alu_txt ALU_TXT [-t THREAD]

Calculate AEI from bam file, both in read level and Alu level

optional arguments:
  -h, --help            show this help message and exit
  -b BAM_FILE, --bam_file BAM_FILE
                        input bam file, with cs tags, sorted and indexed
  --strand_file STRAND_FILE
                        corrrect strand file for reads, generated by the L-GIREMI
  -c [CHROMOSOMES ...], --chromosomes [CHROMOSOMES ...]
                        chromosomes to be analyzed
  -o OUTPUT_PREFIX, --output_prefix OUTPUT_PREFIX
                        prefix of output file
  --genome_fasta GENOME_FASTA
                        path of genome fasta file
  --snp_bcf SNP_BCF     path of dbSNP bcf file
  --alu_txt ALU_TXT     path of bed file of simple repeats [chromosom, start, end, name] (0 based)
  -t THREAD, --thread THREAD
                        cores to be used
```

### `get_read_site`

`get_read_site` output a table that has the read name, genomic
position, and the nucleotide of the position on the read. Running the
script requires a file that provides the genomic position of sites
that are interested.

```{bash}
usage: get_read_site [-h] -s SITE_FILE -b BAM_FILE [-c [CHROMOSOMES ...]] [-o OUTPUT_PREFIX]
                        [--mapq_threshold MAPQ_THRESHOLD] [-t THREAD]

Get read and site information from bam file

optional arguments:
  -h, --help            show this help message and exit
  -s SITE_FILE, --site_file SITE_FILE
                        input site file: [chromosome, pos], without header, separated by tab
  -b BAM_FILE, --bam_file BAM_FILE
                        input bam file, with cs tags, sorted and indexed
  -c [CHROMOSOMES ...], --chromosomes [CHROMOSOMES ...]
                        chromosomes to be analyzed
  -o OUTPUT_PREFIX, --output_prefix OUTPUT_PREFIX
                        prefix of output file
  --mapq_threshold MAPQ_THRESHOLD
                        Min MAPQ to be considered in bam file (default: 20)
  -t THREAD, --thread THREAD
                        cores to be used
```

### `get_read_mismatch`

`get_read_mismatch` is similar to `get_read_site`, except that
it only outputs the read name and the mismatched sites, without matched
sites. `get_read_mismatch` also requires a file that provides the
genomic position of sites that are interested.

```{bash}
usage: get_read_mismatch [-h] -s SITE_FILE -b BAM_FILE [-c [CHROMOSOMES ...]] [-o OUTPUT_PREFIX]
                            [--mapq_threshold MAPQ_THRESHOLD] [-t THREAD]

Get read and mismatch site information from bam file

optional arguments:
  -h, --help            show this help message and exit
  -s SITE_FILE, --site_file SITE_FILE
                        input site file: [chromosome, pos], without header, separated by tab
  -b BAM_FILE, --bam_file BAM_FILE
                        input bam file, with cs tags, sorted and indexed
  -c [CHROMOSOMES ...], --chromosomes [CHROMOSOMES ...]
                        chromosomes to be analyzed
  -o OUTPUT_PREFIX, --output_prefix OUTPUT_PREFIX
                        prefix of output file
  --mapq_threshold MAPQ_THRESHOLD
                        Min MAPQ to be considered in bam file (default: 20)
  -t THREAD, --thread THREAD
                        cores to be used
```

### `get_read_splice`

`get_read_splice` outputs read names and the splicing sites
detected by minimap2 in the reads.

```{bash}
usage: get_read_splice [-h] -b BAM_FILE [-c [CHROMOSOMES ...]] [-o OUTPUT_PREFIX]
                          [--mapq_threshold MAPQ_THRESHOLD] [-t THREAD]

Get read and splice sites from bam file

optional arguments:
  -h, --help            show this help message and exit
  -b BAM_FILE, --bam_file BAM_FILE
                        input bam file, with cs tags, sorted and indexed
  -c [CHROMOSOMES ...], --chromosomes [CHROMOSOMES ...]
                        chromosomes to be analyzed
  -o OUTPUT_PREFIX, --output_prefix OUTPUT_PREFIX
                        prefix of output file
  --mapq_threshold MAPQ_THRESHOLD
                        Min MAPQ to be considered in bam file (default: 20)
  -t THREAD, --thread THREAD
                        cores to be used
```

### `calculate_site_splice_mi`

`calculate_site_splice_mi` can calculate the mutual information for
mismatched sites (RNA editing sites or SNPs) and splicing sites.

```{bash}
usage: calculate_site_splice_mi [-h] [-m READ_SITE] [-s READ_SPLICE] [-o OUTPUT_PREFIX]

calculate the mutual information of mismatch site and splice site pairs

optional arguments:
  -h, --help            show this help message and exit
  -m READ_SITE, --read_site READ_SITE
                        read-site file: [read_name, chromosome, pos, seq]
  -s READ_SPLICE, --read_splice READ_SPLICE
                        corrected read-splice file: [read_name, chromosome, pos, type, corrected_pos,
                        annotation]
  -o OUTPUT_PREFIX, --output_prefix OUTPUT_PREFIX
                        prefix of output file
```

### `split_bam_by_site`

`split_bam_by_site` can fetch reads that cover one genomic location
and save the reads into separate SAM files by the location, which can
be used in IGV ploting.

```{bash}
usage: split_bam_by_site [-h] -b BAM_FILE [-o OUTPUT_PREFIX] [-c CHROMOSOME] [-p POS]

Save reads that cover one position into separated SAM files by the genomic location

optional arguments:
  -h, --help            show this help message and exit
  -b BAM_FILE, --bam_file BAM_FILE
                        input bam file, with cs tags, sorted and indexed
  -o OUTPUT_PREFIX, --output_prefix OUTPUT_PREFIX
                        prefix of output file
  -c CHROMOSOME, --chromosome CHROMOSOME
  -p POS, --pos POS
```

***

[1]: <http://www.nature.com/articles/s41592-019-0610-9> "Roth, S.H., Levanon, E.Y., and Eisenberg, E. (2019). Genome-wide quantification of ADAR adenosine-to-inosine RNA editing activity. Nat Methods 16, 1131–1138."
[2]: <https://academic.oup.com/bioinformatics/article/35/2/340/5038460> "Wyman, D., and Mortazavi, A. (2019). TranscriptClean: variant-aware correction of indels, mismatches and splice junctions in long-read transcripts. Bioinformatics 35, 340–342."
