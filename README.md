# L-GIREMI

L-GIREMI (Long-read Genome-independent Identification of RNA Editing
by Mutual Information) is a software for analysis of RNA editing sites
from long-read RNA-seq data which adopts similar logic with the GIREMI
software [(Zhang and Xiao, 2015)][1].

## Requirement

The L-GIREMI software was developed with python3 on Linux, which
demands several python packages.

* python3.5+: with sys, argparse, re, functools, collections,
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

It's advised to install the software into a virtual environment.

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
python setup.py install
```

#### From pip

```{bash}
pip install l-giremi
```

### Reference data

Be aware to use reference files with same assembly version and the
same chromosome name pattern in all the process.

#### Genome fasta file

The L-GIREMI requires a genome fasta file as input. For human genome,
the fasta file can be obtained from
[UCSC](https://hgdownload.cse.ucsc.edu/goldenpath/hg38/chromosomes/)
or [NCBI](https://www.ncbi.nlm.nih.gov/genome/guide/human/). The fasta
file should include sequences from all the chromosomes.

#### SNP vcf file

The L-GIREMI will use SNP vcf file to get putative heterozygous
SNPs. dbSNP vcf are a possible reference. And it's even better to use
sample specific SNP vcf files. The vcf files should be converted into
bcf format, sorted, and indexed.

The following is a example to convert the dbSNP chromosomes into UCSC chromosomes.

```{bash}
wget "https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.38.gz"

wget "https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.38.gz.md5"

md5sum -c GCF_000001405.38.gz.md5

wget "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.28_GRCh38.p13/GCA_000001405.28_GRCh38.p13_assembly_report.txt"

egrep -v "^#" GCA_000001405.28_GRCh38.p13_assembly_report.txt | cut -f 7,10 | tr "\t" " " > GRCh38-to-hg38.map

bcftools annotate --rename-chrs GRCh38-to-hg38.map GCF_000001405.38.gz -Ob -o dbsnp.38.hg38.bcf

bcftools index dbsnp.38.hg38.bcf
```

#### Gene annotation gtf

A gene annotation gtf is required for the L-GIREMI, which can be
obtained from [gencode](https://www.gencodegenes.org/human/). The gtf
file should be sorted, zipped, and indexed.

The following are codes for gtf file preparation:
```{bash}
(grep ^"#" gencode.v37.annotation.gtf; \
    grep -v ^"#" gencode.v37.annotation.gtf | \
    sort -k1,1 -k4,4n) | \
    bgzip > gencode.v37.annotation.sorted.gtf.gz

tabix -p gff gencode.v37.annotation.sorted.gtf.gz
```

#### Simple repeat region table

The simple repeat region table is used by L-GIREMI to filter sites. It
can be obtained from [UCSC table
browser](http://genome.ucsc.edu/cgi-bin/hgTables). The table format
should be: chromosome, start, end. No header needed for the file.

## Usage

```{bash}
usage: l-giremi [-h] -b BAM_FILE [-c [CHROMOSOMES ...]] [-o OUTPUT_PREFIX] [-t THREAD] --genome_fasta GENOME_FASTA
                  --snp_bcf SNP_BCF --repeat_txt REPEAT_TXT --annotation_gtf ANNOTATION_GTF
                  [--mapq_threshold MAPQ_THRESHOLD] [--min_allele_count MIN_ALLELE_COUNT]
                  [--gene_padding GENE_PADDING] [--exon_padding EXON_PADDING] [--min_rc_cov MIN_RC_COV]
                  [--homopoly_length HOMOPOLY_LENGTH] [--min_AB MIN_AB] [--min_AC MIN_AC]
                  [--min_het_snp_ratio MIN_HET_SNP_RATIO] [--max_het_snp_ratio MAX_HET_SNP_RATIO]
                  [--mi_min_common_read MI_MIN_COMMON_READ] [--mi_min_read MI_MIN_READ] [--mip_threshold MIP_THRESHOLD]

L-GIREMI (Long-read Genome-independent Identification of RNA Editing by Mutual Information)

optional arguments:
  -h, --help            show this help message and exit
  -b BAM_FILE, --bam_file BAM_FILE
                        input bam file, sorted and indexed
  -c [CHROMOSOMES ...], --chromosomes [CHROMOSOMES ...]
                        chromosomes to be analyzed
  -o OUTPUT_PREFIX, --output_prefix OUTPUT_PREFIX
                        prefix of output file
  -t THREAD, --thread THREAD
                        cores to be used
  --genome_fasta GENOME_FASTA
                        path of genome fasta file
  --snp_bcf SNP_BCF     path of dbSNP bcf file
  --repeat_txt REPEAT_TXT
                        path of txt file of simple repeats [chromosom, start, end] (0 based)
  --annotation_gtf ANNOTATION_GTF
                        gtf (gz and tabix indexed) file of genome annotation (gencode)
  --mapq_threshold MAPQ_THRESHOLD
                        Min MAPQ to be considered in bam file (default: 20)
  --min_allele_count MIN_ALLELE_COUNT
                        Min allele read count (default: 2)
  --drop_non_spliced_read DROP_NON_SPLICED_READ
                        Drop non spliced reads (default: True)
  --min_dist_from_splice MIN_DIST_FROM_SPLICE
                        Drop sites within the distance from splice junctions (default: 4)
  --gene_padding GENE_PADDING
                        expand the range when searching gene gtf (default: 500)
  --exon_padding EXON_PADDING
                        expand the range when searching exon gtf (default: 10)
  --min_rc_cov MIN_RC_COV
                        min coverage of read cluster to be considered (default: 2)
  --homopoly_length HOMOPOLY_LENGTH
                        left and right sequence length to be searched for the homopoly around sites (default: 5)
  --min_AB MIN_AB       Min mismatch ratio to be considered (default: 0.1)
  --min_AC MIN_AC       Min mismatch read count to be considered (default: 3)
  --min_het_snp_ratio MIN_HET_SNP_RATIO
                        Min ratio to be considered as heterogenous SNPs (default: 0.35)
  --max_het_snp_ratio MAX_HET_SNP_RATIO
                        Max ratio to be considered as heterogenous SNPs (default: 0.65)
  --mi_min_common_read MI_MIN_COMMON_READ
                        Min common read for site pairs to calculate MI (default: 6)
  --mi_min_read MI_MIN_READ
                        Min read for a variant of a site in a site pair to calculate MI (default: 1)
  --mip_threshold MIP_THRESHOLD
                        MI p value threshold to be used to separate RNA editing sites (default: 0.05)
```

## Analysis process

Running L-GIREMI requires some annotation files:
* reference FASTA file: for example hg38 fasta.
* SNP VCF file: dbSNP VCF file (the chromosome names should be agreed
  with the reference FASTA file), or known SNP VCF file for teh
  sample.
* genome annotation GTF file: with exon annotation and gene
  annotation, for example GENCODE gtf file.

Firstly, map the long-read RNA-seq fastq file using
[minimap2](https://github.com/lh3/minimap2) with cs tag information
provided, which is required by `l-giremi.py`.

```{bash}
minimap2 -t 8 -ax splice -uf \
         --secondary=no -N 5 --cs \
         -o SAM_FILE GENOME_FILE FASTQ_FILE
```

Then, remove unmapped and non-unique mapped reads from SAM file, and
sort SAM into an indexed BAM file.

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

## Output

`l-giremi` generates several result files.

* corrected read strand files: the files are saved my chromosome, and
  are stored as `$OUTPREFIX.$CHROMOSOME.strand`.
  columns:
  1. read_name: the read name.
  2. seq_strand: original read mapping strand.
  3. read_strand: corrected read strand.
* site position files: stored as `$OUTPREFIX.site`.
  columns:
  1. chromosome: chromosome name.
  2. pos: position on the chromosome, 0-based.
  3. ref: sequence of reference genome, strandless.
  4. snp: SNP mark, 0 for not found overlapping with input SNP
     annotations, 1 for overlapping with input SNP annotations.
* mutual information files: stored as `$OUTPREFIX.mi`.
  columns:
  1. type: het_snp, for mismatch sites that overlapping with SNP
     annotations and with mismatch ratio satisfying the
     parameters. mimatch for other types.
  2. chromosome: chromosome name.
  3. pos: position on the chromosome, 0-based.
  4. strand: strand of the mismatch sites.
  5. change_type: the mismatch type, [ref]>[alt].
  6. mi: site mutual information.
  7. n: n sites pairs for the calculation of site mutual information.
  8. pairs: the names of the paired sites.
  9. jakarta: Jakarta index for each site pairs.
  10. mis: mutual information for each site pairs.
  11. mi_cov: mutual information read coverage.
  12. mip: emperical p value of the MI
* mismatch score file: stored as `$OUTPREFIX.score`.
  columns:
  1. type: het_snp, for mismatch sites that overlapping with SNP
     annotations and with mismatch ratio satisfying the
     parameters. mimatch for other types.
  2. chromosome: chromosome name.
  3. pos: position on the chromosome, 0-based.
  4. strand: strand of the mismatch sites.
  5. change_type: the mismatch type, [ref]>[alt].
  6. read_count: the read count for the mismatch.
  7. depth: the total read count for the site.
  8. ratio: the mismatch ratio.
  9. up_seq: a 5' nucleotide ahead of the mismatch sites.
  10. down_seq: a 3' nucleotide after the mismatch sites.
  11. score: RNA editing score by the GLM model.

## Useful scripts

There are also some useful scripts for download anslysis provided with the main l-giremi script.

### `get_aei`

`get_aei.py` is used for the calculation of Alu Editing Index (AEI)
from long-read RNA-seq data [(Roth et al., 2019)][2]. Calculation of
AEI neads the read strand files (output of l-giremi.py), the genome
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

`get_read_mismatch` is similar to `get_read_site.py`, except that
it only output the read name and the mismatched sites, without matched
sites. `get_read_mismatch.py` also requires a file that provides the
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

### `correct_splice_site`

Since the raw splicing sites detected in the long-read RNA-seq data
have many errors, the splicing sites should be corrected [(Wyman and
Mortazavi, 2019)][3]. `correct_splice_site` adopted similar
correction strategies with TranscriptClean. The correction needs a GTF
annotation file.

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
and save the reads into separated SAM files by the location, which can
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

[1]: <http://www.nature.com/articles/nmeth.3314> "Zhang, Q., and Xiao, X. (2015). Genome sequence–independent identification of RNA editing sites. Nat Methods 12, 347–350."
[2]: <http://www.nature.com/articles/s41592-019-0610-9> "Roth, S.H., Levanon, E.Y., and Eisenberg, E. (2019). Genome-wide quantification of ADAR adenosine-to-inosine RNA editing activity. Nat Methods 16, 1131–1138."
[3]: <https://academic.oup.com/bioinformatics/article/35/2/340/5038460> "Wyman, D., and Mortazavi, A. (2019). TranscriptClean: variant-aware correction of indels, mismatches and splice junctions in long-read transcripts. Bioinformatics 35, 340–342."
