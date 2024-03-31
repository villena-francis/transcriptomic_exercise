# Practical assigment for "Transcriptomics" 2023-2024

# Purpose
This repository serves as a submission for the ***Final transcriptomics exercise*** supervised by [@mj-jimenez](https://github.com/mj-jimenez), [@SGMartin](https://github.com/SGMartin) and Jaime Martínez de Villarreal, as part of the Bioinformatics Master's program at ENS-ISCIII.

The original assigment prompt is aviable in this [repository](https://github.com/bioinfo-lessons/transcriptomic-final-exercise)

# Contents
[SECTION 1](#section-1)

[1. Quality control with FastQC, cross-contamination analysis with FastQScreen and read pre-processing with cutadapt](#1-quality-control-with-fastqc-cross-contamination-analysis-with-fastqscreen-and-read-pre-processing-with-cutadapt)

[2. Align samples to chromosome 21 reference with HISAT2](#2-align-samples-to-chromosome-21-reference-with-hisat2)

[3. Aligment statistics and quantify expression using GTF file](#3-aligment-statistics-and-quantify-expression-using-gtf-file)

[SECTION 2](#section-2)

[4. Differential gene expression analysis](#4-differential-gene-expression-analysis)

[5. Gene set enrichment analysis](#5-gene-set-enrichment-analysis)

# Section 1
To address the questions posed in the first section, instructors provided two samples (SRR479052 and SRR479054), consisting of four FASTQ files aviable in [`part_1/data_raw`](part_1/data_raw). They also supplied a FASTA file containing the GRCh38 assembly of human chromosome 21 and a GTF file with corresponding gene annotations (GRCh38.ensembl.109).

The specifications of the mamba environment containing the necessary programs are provided in [`env/part_1.yml`](env/part_1.yml). To create and activate the environment, use the following commands:
```
mamba env create -f env/part_1.yml

mamba activate part_1
```
The [`scripts/section1.sh`](scripts/section1.sh) has been created to facilitate the completion of tasks in this section. The script's configuration and functionality will be explained in detail while addressing the corresponding questions. To execute the pipeline, navigate to the main folder of this repo in your terminal and run the following command:
```
bash scripts/part_1.sh
```
> Running FastQ Screen for the first time requires considerable time due to the need to download large reference files and set the .conf configuration file to indicate the path of the aligner that will be used. For this reason, FastQ Screen is disabled by default in the pipeline. If you want to run it, simply remove the comments (hash symbols #) from the corresponding lines to allow its execution.

## 1. Quality control with FastQC, cross-contamination analysis with FastQScreen and read pre-processing with cutadapt

**FastQC reports** for each fastq file are available in `part_1/quality_controls/fastqc_reports`. Here you can review FastQC reports obtained:

- [SRR479052.chr21_1_fastqc.html](https://htmlpreview.github.io/?https://github.com/villena-francis/transcriptomic_exercise/blob/main/part_1/quality_controls/fastqc_reports/SRR479052.chr21_1_fastqc.html)
- [SRR479052.chr21_2_fastqc.html](https://htmlpreview.github.io/?https://github.com/villena-francis/transcriptomic_exercise/blob/main/part_1/quality_controls/fastqc_reports/SRR479052.chr21_2_fastqc.html)
- [SRR479054.chr21_1_fastqc.html](https://htmlpreview.github.io/?https://github.com/villena-francis/transcriptomic_exercise/blob/main/part_1/quality_controls/fastqc_reports/SRR479054.chr21_1_fastqc.html)
- [SRR479054.chr21_2_fastqc.html](https://htmlpreview.github.io/?https://github.com/villena-francis/transcriptomic_exercise/blob/main/part_1/quality_controls/fastqc_reports/SRR479054.chr21_2_fastqc.html)
  
According to those reports, there are several widespread quality issues with the samples:

- Per base sequence quality. Quality scores drop below 20 in the last 20 positions of the reads. 
- Per base sequence content. nucleotide composition bias along the first 15 positions.
- Overrepresented sequences. In the files corresponding to SRR479052 sample, a sequence matching 18S rRNA is overrepresented.
- Adapter Content. around 10% of the reads still contain Illumina adapter sequences.

**FastQ Screen reports** for each fastq file are available in `part_1/quality_controls/fastq_screen_reports`. Here you can review FastQ Screen reports obtained:
- [SRR479052.chr21_1_screen.html](https://htmlpreview.github.io/?https://github.com/villena-francis/transcriptomic_exercise/blob/main/part_1/quality_controls/fastq_screen_reports/SRR479052.chr21_1_screen.html)
- [SRR479052.chr21_2_screen.html](https://htmlpreview.github.io/?https://github.com/villena-francis/transcriptomic_exercise/blob/main/part_1/quality_controls/fastq_screen_reports/SRR479052.chr21_2_screen.html)
- [SRR479054.chr21_1_screen.html](https://htmlpreview.github.io/?https://github.com/villena-francis/transcriptomic_exercise/blob/main/part_1/quality_controls/fastq_screen_reports/SRR479054.chr21_1_screen.html)
- [SRR479054.chr21_2_screen.html](https://htmlpreview.github.io/?https://github.com/villena-francis/transcriptomic_exercise/blob/main/part_1/quality_controls/fastq_screen_reports/SRR479054.chr21_2_screen.html)

According to these reports, there is no evidence of cross-contamination, although like with FastQC, the existence of sequences related to sequencing adapters is observed. 

Based on the quality control results, it was decided to add a **pre-processing step** with `cutadapt` before aligning the samples. From the [tool's guide](https://cutadapt.readthedocs.io/en/v1.8/guide.html), the following parameters were set:

- `-a AGATCGGAAGAGC -A AGATCGGAAGAGC`: specifies the 3' adapter sequences to trim from the 3' ends of reads 1 and 2, respectively.
 
- `-q 20`: trims low-quality bases from the 3' ends with a quality value below 20.
  
- `-m 50`: Discards reads that have a length less than 50 bases after trimming.
   
- `-u 15 -U 15`: Removes the first 15 bases from reads 1 and 2, respectively.
 
- `-g CTTTTACTTCCTCTAGATAGTCAAGTTCGACCGTCTTCTCAGCGCTCCGC`: removes the 18S rRNA sequence from the reads.
   
- `-o "$trimmed_dir/${sid}.chr21_1_trimmed.fastq"`: output path for trimmed reads 1.

- `-p "$trimmed_dir/${sid}.chr21_2_trimmed.fastq"`: output path for trimmed reads 2.

- `$fastq_dir/${sid}.chr21_1.fastq $fastq_dir/${sid}.chr21_2.fastq`: input paths for FASTQ files of reads 1 and 2, respectively.

- `&> $trimmed_dir/cutadapt_stats.log`: save the cutadapt execution statistics in a log file.
  
Trimmed files are aviable in [`part_1/cutadapt`](part_1/cutadapt).

To ensure the pre-processing steps were successful, FastQC analysis was performed on the trimmed FASTQ files. These reports are stored in `part_1/quality_control/fastqc_reports` and are labeled with "trimmed" to distinguish them from the reports of the raw data. Here you can review FastQC reports obtained from the trimmed FASTQ files:

- [SRR479052.chr21_1_trimmed_fastqc.html](https://htmlpreview.github.io/?https://github.com/villena-francis/transcriptomic_exercise/blob/main/part_1/quality_controls/fastqc_reports/SRR479052.chr21_1_trimmed_fastqc.html)
- [SRR479052.chr21_2_trimmed_fastqc.html](https://htmlpreview.github.io/?https://github.com/villena-francis/transcriptomic_exercise/blob/main/part_1/quality_controls/fastqc_reports/SRR479052.chr21_2_trimmed_fastqc.html)
- [SRR479054.chr21_1_trimmed_fastqc.html](https://htmlpreview.github.io/?https://github.com/villena-francis/transcriptomic_exercise/blob/main/part_1/quality_controls/fastqc_reports/SRR479054.chr21_1_trimmed_fastqc.html)
- [SRR479054.chr21_2_trimmed_fastqc.html](https://htmlpreview.github.io/?https://github.com/villena-francis/transcriptomic_exercise/blob/main/part_1/quality_controls/fastqc_reports/SRR479054.chr21_2_trimmed_fastqc.html)

The outcomes of the reports were favorable across all parameters, except for a warning in Sequence Length Distribution. This warning about the increased variability in sequence lengths was expected after cutadapt pre-processing.

## 2. Align samples to chromosome 21 reference with HISAT2 

To **align sequencing reads with HISAT2**, it is first necessary to build an index of the reference genome. This is done using the `hisat2-build` command and the following arguments:

- `--seed 123`: sets a seed for the random number generator, allowing for reproducible results.
- `-p 10`: indicates that 10 threads should be used for index construction.
- `genomes/*chromosome.21.fa`: specifies the FASTA file containing the sequence of chromosome 21.
- `genomes/hg38_chr21`: specifies the prefix for the generated index files.
- `&> genomes/hisat2_index.log`: save the hisat2-build execution statistics in a log file.

Here you can review the [HISAT2 index statistics](part_1/genomes/hisat2_index.log).

Once the reference index has been built, sequencing reads can be aligned using the `hisat2` command with the following arguments:

- `--new-summary`: Indicates that a new summary file should be created.
- `--summary-file ${aligments_dir}/hisat2_${sid}.stats`: specifies the path and name of the generated stats file.
- `--seed 123`: sets a seed for the random number generator, allowing for reproducible results.
- `--phred33`: indicates that the quality scores in the FASTQ files are encoded in Phred+33 format.
- `-p 10`: indicates that 10 threads should be used for alignment, speeding up the process.
- `-k 1`: indicates that only the best alignment for each read should be reported.
- `-x genomes/hg38_chr21`: specifies the prefix of the index files previously generated with hisat2-build.
- `-1 ${trimmed_dir}/${sid}.chr21_1_trimmed.fastq`: specifies the FASTQ file containing the trimmed forward reads.
- `-2 ${trimmed_dir}/${sid}.chr21_2_trimmed.fastq`: specifies the FASTQ file containing the trimmed reverse reads.
- `-S ${aligments_dir}/${sid}_trimmed.sam`: specifies the path and name of the output SAM file that will contain the alignments.

After aligning the reads with HISAT2, the resulting SAM files were processed using the **samtools toolkit** to convert to BAM format, sort the alignments and index the BAM files.

## 3. Aligment statistics and quantify expression using GTF file

The statistics of the alignments were generated with the execution of the `hisat2` command and are available in `part_1/hisat2` as `.stats` files. Here you can review the HISAT2 aligment statistics by sample:
- [hisat2_SRR479052.stats](part_1/hisat2/hisat2_SRR479052.stats)
- [hisat2_SRR479054.stats](part_1/hisat2/hisat2_SRR479054.stats)

Finally, **gene expression levels** were quantified using the `htseq-count` tool with a GTF annotation file using the following arguments:

- `--format=bam`: indicates that the input file is in BAM format.
- `--stranded=reverse`: indicates that the reads come from a reverse strand-specific RNA-seq library.
- `--mode=intersection-nonempty`: specifies the counting mode, where a read is counted if it overlaps with at least one base of an exon and does not extend beyond its boundaries.
- `--minaqual=10`: sets a minimum alignment quality threshold of 10 for a read to be counted.
- `--type=exon`: specifies that reads overlapping with exons should be counted.
- `--idattr=gene_id --additional-attr=gene_name`: indicates that the "gene_id" attribute should be used as the feature ID and that the "gene_name" attribute should also be included in the output.
- `${aligments_dir}/${sid}_trimmed_sorted.bam`: specifies the input BAM file with the processed alignments.
- `genomes/Homo_sapiens.GRCh38.109.chr21.gtf`: specifies the GTF annotation file containing the genomic feature coordinates.
- `> ${quant_dir}/${sid}.htseq`: redirects the output to the specified file, which will contain the read counts per gene.

The obtained read counts are available in `part_1/htseq` as `.htseq` files:
- [SRR479052.htseq](part_1/htseq/SRR479052.htseq)
- [SRR479054.htseq](part_1/htseq/SRR479054.htseq)