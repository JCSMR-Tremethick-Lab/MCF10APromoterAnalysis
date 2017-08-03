#! /bin/bash
# Arg1 - Arrowhead domainlist
# Arg2 - chromosome size --- Make sure it is sorted and chrY removed


filename="${1%.*}"
tail $1 | sort -k1,1 -k2,2n - | mergeBed -i - | awk '{print "chr"$1 "\t" $2 "\t" $3 "\tdomain"}' - > ${filename}.merged.bed
sort -k1,1 -k2,2n $2 > genome
complementBed -i ${filename}.merged.bed -g genome | cat - ${filename}.merged.bed | sort -k1,1 -k2,2n - | awk '{if($4=="domain") print $0; else print $1 "\t" $2 "\t" $3 "\tgap"}'  > ${filename}.domains
