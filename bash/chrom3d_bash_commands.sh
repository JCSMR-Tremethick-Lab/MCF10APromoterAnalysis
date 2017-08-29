#! /bin/bash

cat intra_chr_bedpe/*.bedpe > intra_chr_bedpe/shZ-rep2.bedpe

curl -s "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz" | gunzip -c | grep acen |\
~/miniconda3/envs/chrom3d/bin/pairToBed -a intra_chr_bedpe/shZ-rep2.bedpe -b stdin -type  neither >\
intra_chr_bedpe/shZ-rep2.domain.RAW.no_cen.bedpe

NCHG -m 40000 -p intra_chr_bedpe/shZ-rep2.domain.RAW.no_cen.bedpe > intra_chr_bedpe/shZ-rep2.domain.RAW.no_cen.NCHG.out

# should yield ~ 3000-8000 significant interactions
Rscript ~/Development/JCSMR-Tremethick-Lab/Breast/R/NCHG_fdr_oddratio_calc.R intra_chr_bedpe/shZ-rep2.domain.RAW.no_cen.NCHG.out \
intra_chr_bedpe/shZ-rep2.domain.RAW.no_cen.NCHG.sig 2 BH 0.01

~/Development/JCSMR-Tremethick-Lab/Breast/bash/make_gtrack.sh intra_chr_bedpe/shZ-rep2.domain.RAW.no_cen.NCHG.sig 40k_intra/shZ-rep2.domains shZ-rep2_intra_chromosome.gtrack

~/Development/JCSMR-Tremethick-Lab/Breast/bash/make_gtrack_incl_LADs.sh shZ-rep2_intra_chromosome.gtrack ~/Data/chrom3d_pipeline_data/GSM1313399_HSF_AD04_LMNA_rep1.bed rep2_intra_chromosome_w_LADs.gtrack

mkdir interchr_bedpe

# this is implemented in snakemake
~/miniconda3/envs/chrom3d/bin/python ~/Development/JCSMR-Tremethick-Lab/Breast/Python/make_interchr_NCHG_input.py 1M_inter_chr_RAWobserved/shZ-rep2.6.19.1000k.coo \
~/Data/References/Annotations/Homo_sapiens/hg19/UCSC/unmappable_blacklist.bed \
~/Data/References/Annotations/Homo_sapiens/hg19/UCSC/hg19.chrom.sizes.sorted chr6 chr19 > inter_chr_bedpe/shZ-rep2.6.19.1000k.bedpe

cat inter_chr_bedpe/*.bedpe >> inter_chr_bedpe/shZ-rep2.1MB.bedpe

mv inter_chr_bedpe/shZ-rep2.bedpe inter_chr_bedpe/shZ-rep2.1MB.bedpe
NCHG -i -p inter_chr_bedpe/shZ-rep2.bedpe > inter_chr_bedpe/shZ-rep2.NCHG.out
mv inter_chr_bedpe/shZ-rep2.NCHG.out inter_chr_bedpe/shZ-rep2.NCHG.1MB.out

Rscript ~/Development/JCSMR-Tremethick-Lab/Breast/R/NCHG_fdr_oddratio_calc.R inter_chr_bedpe/shZ-rep2.NCHG.1MB.out inter_chr_bedpe/shZ-rep2.NCHG.1MB.sig 2 BH 0.01

source activate chrom3d
~/Development/JCSMR-Tremethick-Lab/Breast/bash/add_inter_chrom_beads.sh shZ-rep2_intra_chromosome_w_LADs.gtrack \
inter_chr_bedpe/shZ-rep2.NCHG.1MB.sig shZ-rep2_inter_intra_chr_w_LADs.gtrack

~/miniconda3/envs/chrom3d/bin/python ~/Development/JCSMR-Tremethick-Lab/Breast/Python/make_diploid_gtrack.py shZ-rep2_inter_intra_chr_w_LADs.gtrack > shZ-rep2_inter_intra_chr_w_LADs.diploid.gtrack

Chrom3D -y 0.15 -r 5.0 -n 2000000 -o shZ-rep2_inter_intra_chr_w_LADs.diploid.cmm shZ-rep2_inter_intra_chr_w_LADs.diploid.gtrack

# Make HiC data juicebox compatible
samps="shZ-rep2_S2.hdf5"
for s in $samps 
do
nohup python2.7 makeJuice.py $s > $s.txt 2> $s.txt.err &
wait
done
