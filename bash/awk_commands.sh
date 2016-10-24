find mapped-hg38 -name "MCF_Cancer*" | awk -v OFS="\t" '{cell="MCF10A"; replicate="R1"; genome="hg38"; enzyme="HindIII"; print $1"\t"cell"\t"replicate"\t" genome"\t"enzyme}'

find mapped-hg38 -name "MCF_WT*" | awk -v FS="_" -v OFS="\t" '{if ($3==1) replicate="R1"; else if ($3==2) replicate="R2"; cell="MCF10A_WT"; genome="hg38"; enzyme="HindIII"; print $0"\t"cell"\t"replicate"\t" genome"\t"enzyme}' > runs.tsv

find mapped-hg38 -name "MCF_WT*" | awk -v IFS="_" -v OFS="\t" '{if ($3==1) replicate="R1"; print replicate}'
