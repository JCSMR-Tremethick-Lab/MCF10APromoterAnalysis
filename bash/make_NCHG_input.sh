#! /bin/bash
awk -v chr_name=$3 '{printf("%s\t%s\t%s\n", chr_name, $1,$1+1)}' $2 | intersectBed -wao -a stdin -b $1 | cut -f 4,5,6 > left.tmp

awk -v chr_name=$3 '{printf("%s\t%s\t%s\t%s\n", chr_name, $2,$2+1,$3)}' $2 | intersectBed -wao -a stdin -b $1 | awk '{printf("%s\t%s\t%s\t%s\n",$5,$6,$7,$4)}' > right.tmp

paste left.tmp right.tmp | awk '{a[$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6] += $7} END{for (i in a) print i"\t"a[i]}' |  awk '$2!=$5' | sort -k 2n,2n

rm left.tmp right.tmp
