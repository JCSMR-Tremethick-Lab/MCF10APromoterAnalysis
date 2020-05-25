#!/bin/bash
eval "$(conda shell.bash hook)"

function SubSample {

## see also: http://crazyhottommy.blogspot.com/2016/05/downsampling-for-bam-files-to-certain.html
FACTOR=$(samtools idxstats $1 | cut -f3 | awk -v COUNT=$2 'BEGIN {total=0} {total += $1} END {print COUNT/total}')

if [[ $FACTOR > 1 ]]
  then 
  echo '[ERROR]: Requested number of reads exceeds total read count in' $1 '-- exiting' && exit 1
fi

sambamba view -s $FACTOR -t 2 -f bam -l 5 $1

}

export -f SubSample

export genome_chip_dir="/home/sebastian/Data/Tremethick/Breast/ChIP-Seq/merged/processed_data/GRCh37_hg19_ensembl75_ERCC/samtools/merge/duplicates_removed"
export capture_chip_dir="/home/sebastian/Data/Tremethick/Breast/PromoterSeqCap/bam"

mkdir -p /home/sebastian/Data/Tremethick/Breast/sub_sampling_analysis

# subsample BAM files to roughly equivalent numbers of reads/library
cd $genome_chip_dir
ls *.bam | grep -v subsampled | parallel -j 2 "SubSample {} 5000000 > {.}_subsampled_5M.bam"
ls *.bam | grep -v subsampled | parallel -j 2 "SubSample {} 10000000 > {.}_subsampled_10M.bam"
ls *.bam | grep -v subsampled | parallel -j 2 "SubSample {} 20000000 > {.}_subsampled_20M.bam"

## Capture-seq data
cd $capture_chip_dir/subsampling 
ls *.bam | grep -v subsampled | parallel -j 2 "SubSample {} 5000000 > {.}_subsampled_5M.bam"
ls *.bam | grep -v subsampled | parallel -j 2 "SubSample {} 10000000 > {.}_subsampled_10M.bam"
ls *.bam | grep -v subsampled | parallel -j 2 "SubSample {} 20000000 > {.}_subsampled_20M.bam"


export canonical_chromosomes="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM"
# have to re-header BAM files to match genome annotations
ls *subsampled_5M* | parallel "samtools reheader header.sam {} "

for i in $(ls *subsampled*); do samtools reheader header.sam $i > temp.bam && mv temp.bam $i; done

# 
conda activate deeptools
bamCoverage


