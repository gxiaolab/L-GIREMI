# GIREMIL

GIREMIL (Genome-independent Identification of RNA Editing by Mutual
Information for Long-read RNA-seq)

## Requirement

The GIREMIL software was developed with python3 on Linux, which demands
several python packages.

* python3.5+: with sys, argparse, re, functools, collections, multiprocessing.
* [scikit-learn](https://scikit-learn.org): 0.20+
* [numpy](https://numpy.org): 1.10+
* [pandas](https://pandas.pydata.org): 1.0+
* [pysam](https://pysam.readthedocs.io): 0.16+

And the analysis with GIREMIL also needs additional software:
* [samtools](http://www.htslib.org)
* [minimap2](https://github.com/lh3/minimap2)

## Usage

```{bassh}
usage: giremil.py [-h] -b BAM_FILE [-c [CHROMOSOMES ...]] [-o OUTPUT_PREFIX] [-t THREAD] --genome_fasta GENOME_FASTA
                  --snp_bcf SNP_BCF --repeat_txt REPEAT_TXT --annotation_gtf ANNOTATION_GTF
                  [--mapq_threshold MAPQ_THRESHOLD] [--min_allele_count MIN_ALLELE_COUNT]
                  [--gene_padding GENE_PADDING] [--exon_padding EXON_PADDING] [--min_rc_cov MIN_RC_COV]
                  [--homopoly_length HOMOPOLY_LENGTH] [--min_AB MIN_AB] [--min_AC MIN_AC]
                  [--min_het_snp_ratio MIN_HET_SNP_RATIO] [--max_het_snp_ratio MAX_HET_SNP_RATIO]
                  [--mi_min_common_read MI_MIN_COMMON_READ] [--mi_min_read MI_MIN_READ] [--mi_threshold MI_THRESHOLD]

GIREMIL (Genome-independent Identification of RNA Editing by Mutual Information for Long-read RNA-seq)

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
  --mi_threshold MI_THRESHOLD
                        MI threshold to be used to separate RNA editing sites (default: 0.2)
```

## Analysis process

Running GIREMIL requires some annotation files:
* reference FASTA file: for example hg38 fasta.
* SNP VCF file: dbSNP VCF file (the chromosome names should be agreed
  with the reference FASTA file), or known SNP VCF file for teh
  sample.
* genome annotation GTF file: with exon annotation and gene
  annotation, for example GENCODE gtf file.

Firstly, map the long-read RNA-seq fastq file using
[minimap2](https://github.com/lh3/minimap2) with cs tag information
provided, which is required by `giremil.py`.

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

Next, run `giremil.py`.

```{bash}
giremil.py \
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

`giremil.py` generates several result files.

* corrected read strand files: the files are saved my chromosome, and
  are stored as `$OUTPREFIX.$CHROMOSOME.strand`.
  columns:
  1. read_name: the read name.
  2. seq_strand: original read mapping strand.
  3. read_strand: corrected read strand.
* site position files: the files are saved by chromosome, and are
  stored as `$OUTPREFIX.$CHROMOSOME.site`.
  columns:
  1. chromosome: chromosome name.
  2. pos: position on the chromosome, 0-based.
  3. ref: sequence of reference genome, strandless.
  4. snp: SNP mark, 0 for not found overlapping with input SNP
     annotations, 1 for overlapping with input SNP annotations.
* mutual information files: the files are saved by chromosome, and are
  stored as `$OUTPREFIX.$CHROMOSOME.mi`.
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
* mismatch ratio files: the files are saved by chromosome, and are
  stored as `$OUTPREFIX.$CHROMOSOME.mismatch_ratio`.
  columns:
  1. type: het_snp, for mismatch sites that overlapping with SNP
     annotations and with mismatch ratio satisfying the
     parameters. mimatch for other types.
  2. chromosome: chromosome name.
  3. pos: position on the chromosome, 0-based.
  4. strand: strand of the mismatch sites.
  5. change_type: the mismatch type, [ref]>[alt].
  6. ratio: the mismatch ratio.
* mismatch score file: one file for the analysis, and is stored as
  `$OUTPREFIX.score`.
  columns:
  1. type: het_snp, for mismatch sites that overlapping with SNP
     annotations and with mismatch ratio satisfying the
     parameters. mimatch for other types.
  2. chromosome: chromosome name.
  3. pos: position on the chromosome, 0-based.
  4. strand: strand of the mismatch sites.
  5. change_type: the mismatch type, [ref]>[alt].
  6. ratio: the mismatch ratio.
  7. up_seq: a 5' nucleotide ahead of the mismatch sites.
  8. down_seq: a 3' nucleotide after the mismatch sites.
  9. score: RNA editing score by the GLM model.
