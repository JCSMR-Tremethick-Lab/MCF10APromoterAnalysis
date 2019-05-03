#!/bin/bash
export genome_chip_dir="/home/sebastian/Data/Tremethick/Breast/ChIP-Seq/merged/processed_data/GRCh37_hg19_ensembl75_ERCC/samtools/merge/duplicates_removed"
export capture_chip_dir="/home/sebastian/Data/Tremethick/Breast/PromoterSeqCap/bam"

mkdir -p /home/sebastian/Data/Tremethick/Breast/sub_sampling_analysis

# subsample BAM files to roughly equivalent numbers of reads/library
cd $genome_chip_dir
samtools view -b -s 42.006 Inp_10A_WT_high.bam 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y MT > Inp_10A_WT_high.subsampled.bam &&\
samtools index Inp_10A_WT_high.subsampled.bam &&\
samtools view -c Inp_10A_WT_high.subsampled.bam #2772794
samtools view -b -s 42.01 H2AZ_10A_high.bam 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y MT > H2AZ_10A_high.subsampled.bam &&\
samtools index H2AZ_10A_high.subsampled.bam &&\
samtools view -c H2AZ_10A_high.subsampled.bam #2777006

## Capture-seq data
cd $capture_chip_dir
samtools sort -@ 4 -m 4G A_Inp_H_r2_R1.fastq_q20.bam -f MCF10A_Input_H_r2_sorted.bam
samtools sort -@ 4 -m 4G A_Inp_H_r1_R1.fastq_q20.bam -f MCF10A_Input_H_r1_sorted.bam
samtools sort -@ 4 -m 4G A_H2AZ_H_r1_R1.fastq_q20.bam -f MCF10A_H2AZ_H_r1_sorted.bam
samtools sort -@ 4 -m 4G A_H2AZ_H_r2_R1.fastq_q20.bam -f MCF10A_H2AZ_H_r2_sorted.bam

export canonical_chromosomes="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM"
# have to re-header BAM files to match genome annotations
samtools view -b -s 42.10 MCF10A_H2AZ_H_r1_sorted.bam $(echo $canonical_chromosomes) > MCF10A_H2AZ_H_r1_sorted.subsampled.bam &&\
samtools index MCF10A_H2AZ_H_r1_sorted.subsampled.bam &&\
samtools view -c MCF10A_H2AZ_H_r1_sorted.subsampled.bam #2757580
samtools reheader header.sam MCF10A_H2AZ_H_r1_sorted.subsampled.bam  > temp.bam && mv temp.bam MCF10A_H2AZ_H_r1_sorted.subsampled.bam

samtools view -b -s 42.10 MCF10A_H2AZ_H_r2_sorted.bam $(echo $canonical_chromosomes) > MCF10A_H2AZ_H_r2_sorted.subsampled.bam &&\
samtools index MCF10A_H2AZ_H_r2_sorted.subsampled.bam &&\
samtools view -c MCF10A_H2AZ_H_r2_sorted.subsampled.bam &&\
samtools reheader header.sam MCF10A_H2AZ_H_r2_sorted.subsampled.bam  > temp.bam && mv temp.bam MCF10A_H2AZ_H_r2_sorted.subsampled.bam

samtools view -b -s 42.15 MCF10A_Input_H_r1_sorted.bam $(echo $canonical_chromosomes) > MCF10A_Input_H_r1_sorted.subsampled.bam &&\
samtools index MCF10A_Input_H_r1_sorted.subsampled.bam &&\
samtools view -c MCF10A_Input_H_r1_sorted.subsampled.bam &&\
samtools reheader header.sam MCF10A_Input_H_r1_sorted.subsampled.bam  > temp.bam && mv temp.bam MCF10A_Input_H_r1_sorted.subsampled.bam

samtools view -b -s 42.20 MCF10A_Input_H_r2_sorted.bam $(echo $canonical_chromosomes) > MCF10A_Input_H_r2_sorted.subsampled.bam &&\
samtools index MCF10A_Input_H_r2_sorted.subsampled.bam &&\
samtools view -c MCF10A_Input_H_r2_sorted.subsampled.bam &&\
samtools reheader header.sam MCF10A_Input_H_r2_sorted.subsampled.bam  > temp.bam && mv temp.bam MCF10A_Input_H_r2_sorted.subsampled.bam

# 
conda activate deeptools
bamCoverage


