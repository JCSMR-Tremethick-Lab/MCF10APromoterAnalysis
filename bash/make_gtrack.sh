#! /bin/sh
~/miniconda3/envs/chrom3d/bin/python ~/Development/JCSMR-Tremethick-Lab/Breast/Python/make_gtrack.py $1 $2 | sort -k1,1 -k2,2n > $3
sed -i  '1s/^/##gtrack version: 1.0\n##track type: linked segments\n###seqid\tstart\tend\tid\tradius\tedges\n/' $3
