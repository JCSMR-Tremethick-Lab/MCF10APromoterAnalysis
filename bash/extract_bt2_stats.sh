#!/bin/bash

((amb=0));
((exact=0));
((total=0));
for i in *global_hg38_ensembl84.log;
  do reads=$(head -1 $i| cut -f 1 -d " ");
    sample=$(pwd| rev | cut -d'/' -f 1 | rev)
    exact=$(head -4 $i|tail -1| cut -f 5 -d " ");
    amb=$(head -5 $i|tail -1| cut -f 5 -d " ");
    ((total+=$reads));
    ((exact+=$exact));
    ((amb+=$amb));
done
echo -e $sample "\t" "total: "$total "\t" "exact: "$exact "\t" "multi: "$amb;
