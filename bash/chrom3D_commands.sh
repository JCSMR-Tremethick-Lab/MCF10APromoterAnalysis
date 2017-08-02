sh ./make_NCHG_input.sh GSE63525_IMR90_Arrowhead_domainlist.chr19.domains \
intra_chr_RAWobserved/chr19_50kb.RAWobserved chr19 > \
intrachr_bedpe/chr19_50kb.domains.RAW.bedpe

sh ./intrachr_NCHG_input_auto.sh hg19.chrom.sizes.sorted

cat intrachr_bedpe/chr*.bedpe > intrachr_bedpe/IMR90_50kb.domain.RAW.bedpe

curl -s "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz" | gunzip -c | grep acen | pairToBed -a \
intrachr_bedpe/IMR90_50kb.domain.RAW.bedpe -b stdin -type  neither > \
intrachr_bedpe/IMR90_50kb.domain.RAW.no_cen.bedpe

NCHG -m 50000 -p intrachr_bedpe/IMR90_50kb.domain.RAW.no_cen.bedpe > \
IMR90_50kb.domain.RAW.no_cen.NCHG.out

Rscript NCHG_fdr_oddratio_calc.R IMR90_50kb.domain.RAW.no_cen.NCHG.out \
IMR90_50kb.domain.RAW.no_cen.NCHG.sig 2 BH 0.01

sh ./make_gtrack.sh IMR90_50kb.domain.RAW.no_cen.NCHG.sig \
GSE63525_IMR90_Arrowhead_domainlist.domains \
IMR90_intra_chromosome.gtrack

sh ./make_gtrack_incl_lad.sh IMR90_intra_chromosome.gtrack \ GSM1313399_HSF_AD04_LMNA_rep1.bed IMR90_intra_chromosome_w_LADs.gtrack

python make_interchr_NCHG_input.py \
inter_chr_RAWobserved/chr14_15_1mb.RAWobserved \
unmappable_blacklist.bed hg19.chrom.sizes.sorted chr14 chr15 > \
interchr_bedpe/chr14_chr15_1mb.bedpe

sh ./interchr_NCHG_input_auto.sh hg19.chrom.sizes.sorted

cat interchr_bedpe/*chr*.bedpe > interchr_bedpe/IMR90_1mb_inter.bedpe

NCHG -i -p interchr_bedpe/IMR90_1mb_inter.bedpe > IMR90_1mb_inter_chr.NCHG.out

Rscript NCHG_fdr_oddratio_calc.R IMR90_1mb_inter_chr.NCHG.out \
IMR90_1mb_inter_chr.NCHG.sig 2 BH 0.01

sh ./add_inter_chrom_beads.sh IMR90_intra_chromosome_w_LADs.gtrack \
IMR90_1mb_inter_chr.NCHG.sig IMR90_inter_intra_chr_w_LADs.gtrack

python make_diploid_gtrack.py IMR90_inter_intra_chr_w_LADs.gtrack > \
IMR90_inter_intra_chr_w_LADs.diploid.gtrack

Chrom3D -y 0.15 -r 5.0 -n 2000000 -o \
IMR90_inter_intra_chr_w_LADs.diploid.cmm \
IMR90_inter_intra_chr_w_LADs.diploid.gtrack
