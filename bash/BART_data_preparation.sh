#!/bin/sh
export liftOver="/home/sebastian/Bioinformatics/UCSC/liftOver"
export chain="/home/sebastian/Data/References/Annotations/Homo_sapiens/hg38/UCSC/hg19ToHg38.over.chain"

cd /home/sebastian/Data/Collaborations/FSU/PromoterSeqCap/SmallFragments

for i in  $(ls *.bed| cut -f 1 -d "."); do echo "$liftOver $i.bed $chain hg38/${i}.bed hg38/${i}.unmapped.bed"; done
