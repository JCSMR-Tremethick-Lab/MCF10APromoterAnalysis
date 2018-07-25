#!/bin/bash

# using sequences under peak
# macs2 call to produce peaks
macs2 callpeak -t TOTALcombined_A_H2AZ_000-125_chrMremoved.bed -f BED --outdir macs2PeakCalling/TOTAL_H2AZ -g 40999507 --nomodel --extsize 147 --trackline --name "WT_total_H2AZ_final_chrMremoved" --bdg --call-summits --mfold 3 60 --qvalue 0.015
meme-chip -oc meme_summitTest -dna -desc "HOCOMOCOv11_full_HUMAN_mono_meme_format.meme" -db ~/Data/References/MEME/motif_databases/HUMAN/ -meme-minw 6 -meme-maxw 12 -meme-nmotifs 6 -centrimo-local -centrimo-maxreg 50 peakCallingTest.fa


# trying different parameters
macs2 callpeak -t TOTALcombined_A_Inp_000-125.bed -f BED --outdir macs2PeakCalling/MCF10A_WT_TOTAL_Input -g 40999507 --trackline --name "WT_total_Input" --bdg --call-summits --mfold 2 80 -q 0.1 --nomodel --extsize 125

macs2 callpeak -t TOTALcombined_A_Inp_000-100.bed -f BED --outdir macs2PeakCalling/MCF10A_WT_TOTAL_Input_100bp -g 40999507 --trackline --name "WT_total_Input_100bp" --bdg --call-summits --mfold 2 80 -q 0.1 --nomodel --extsize 100

# optimising parameters for input Data
cd /home/sebastian/Data/Collaborations/FSU/PromoterSeqCap/SmallFragments

# shifting p-value threshold so that more peaks are detected
macs2 callpeak -t TOTALcombined_A_Inp_000-125.bed -f BED --outdir macs2PeakCalling/MCF10A_WT_TOTAL_Input_125bp -g 40999507 --trackline --name "WT_total_Input_125bp" --bdg --call-summits --nomodel --shift -50 --extsize 100 -p 0.05
cd /home/sebastian/Data/Collaborations/FSU/PromoterSeqCap/SmallFragments/macs2PeakCalling/MCF10A_WT_TOTAL_Input_125bp
# make bigwig files
tail -$(expr $(wc -l WT_total_Input_125bp_treat_pileup.bdg|cut -f 1 -d " ") - 1) WT_total_Input_125bp_treat_pileup.bdg | sort -k1,1 -k2,2n > WT_total_Input_125bp_treat_pileup.sorted.bdg
bedGraphToBigWig WT_total_Input_125bp_treat_pileup.sorted.bdg ~/Data/References/Genomes/Homo_sapiens/GRCh37_hg19_UCSC/chromSizes.txt WT_total_Input_125bp_treat_pileup.sorted.bw

tail -$(expr $(wc -l WT_total_Input_125bp_control_lambda.bdg|cut -f 1 -d " ") - 1) WT_total_Input_125bp_control_lambda.bdg | sort -k1,1 -k2,2n > WT_total_Input_125bp_control_lambda.sorted.bdg
bedGraphToBigWig WT_total_Input_125bp_control_lambda.sorted.bdg ~/Data/References/Genomes/Homo_sapiens/GRCh37_hg19_UCSC/chromSizes.txt WT_total_Input_125bp_control_lambda.sorted.bw

# for Input from TGFb
# first run cutoff-analysis
macs2 callpeak -t TOTALcombined_A_TGFb_H2AZ_000-125.bed -f BED --outdir macs2PeakCalling/MCF10A_TGFb_TOTAL_Input_125bp -g 40999507 --trackline --name "TGFb_total_Input_125bp" --bdg --call-summits --nomodel --shift -50 --extsize 100 --cutoff-analysis
# standard parameters seems to return reasonable number of peaks
# macs2 callpeak -t TOTALcombined_A_TGFb_H2AZ_000-125.bed -f BED --outdir macs2PeakCalling/MCF10A_TGFb_TOTAL_Input_125bp -g 40999507 --trackline --name "TGFb_total_Input_125bp" --bdg --call-summits --nomodel --shift -50 --extsize 100 -p 0.05
cd /home/sebastian/Data/Collaborations/FSU/PromoterSeqCap/SmallFragments/macs2PeakCalling/MCF10A_TGFb_TOTAL_Input_125bp
# make bigwig files
tail -$(expr $(wc -l TGFb_total_Input_125bp_treat_pileup.bdg|cut -f 1 -d " ") - 1) TGFb_total_Input_125bp_treat_pileup.bdg | sort -k1,1 -k2,2n > TGFb_total_Input_125bp_treat_pileup.sorted.bdg
bedGraphToBigWig TGFb_total_Input_125bp_treat_pileup.sorted.bdg ~/Data/References/Genomes/Homo_sapiens/GRCh37_hg19_UCSC/chromSizes.txt TGFb_total_Input_125bp_treat_pileup.sorted.bw

tail -$(expr $(wc -l TGFb_total_Input_125bp_control_lambda.bdg|cut -f 1 -d " ") - 1) TGFb_total_Input_125bp_control_lambda.bdg | sort -k1,1 -k2,2n > TGFb_total_Input_125bp_control_lambda.sorted.bdg
bedGraphToBigWig TGFb_total_Input_125bp_control_lambda.sorted.bdg ~/Data/References/Genomes/Homo_sapiens/GRCh37_hg19_UCSC/chromSizes.txt TGFb_total_Input_125bp_control_lambda.sorted.bw
# re-run with less stringent cutoff
macs2 callpeak -t TOTALcombined_A_TGFb_H2AZ_000-125.bed -f BED --outdir macs2PeakCalling/MCF10A_TGFb_TOTAL_Input_125bp -g 40999507 --name "TGFb_total_Input_125bp" --call-summits --nomodel --extsize 100 -p 0.1



# for Input from shH2AZ
# first run cutoff-analysis
macs2 callpeak -t TOTALcombined_shH2AZ_Inp_000-125.bed -f BED --outdir macs2PeakCalling/MCF10A_shH2AZ_TOTAL_Input_125bp -g 40999507 --trackline --name "shH2AZ_total_Input_125bp" --bdg --call-summits --nomodel --shift -50 --extsize 100 --cutoff-analysis
# changing p-value cut-off to include more peaks
macs2 callpeak -t TOTALcombined_shH2AZ_Inp_000-125.bed -f BED --outdir macs2PeakCalling/MCF10A_shH2AZ_TOTAL_Input_125bp -g 40999507 --trackline --name "shH2AZ_total_Input_125bp" --bdg --call-summits --nomodel --shift -50 --extsize 100 -p 0.1
cd /home/sebastian/Data/Collaborations/FSU/PromoterSeqCap/SmallFragments/macs2PeakCalling/MCF10A_shH2AZ_TOTAL_Input_125bp
# make bigwig files
tail -$(expr $(wc -l shH2AZ_total_Input_125bp_treat_pileup.bdg|cut -f 1 -d " ") - 1) shH2AZ_total_Input_125bp_treat_pileup.bdg | sort -k1,1 -k2,2n > shH2AZ_total_Input_125bp_treat_pileup.sorted.bdg
bedGraphToBigWig shH2AZ_total_Input_125bp_treat_pileup.sorted.bdg ~/Data/References/Genomes/Homo_sapiens/GRCh37_hg19_UCSC/chromSizes.txt shH2AZ_total_Input_125bp_treat_pileup.sorted.bw

tail -$(expr $(wc -l shH2AZ_total_Input_125bp_control_lambda.bdg|cut -f 1 -d " ") - 1) shH2AZ_total_Input_125bp_control_lambda.bdg | sort -k1,1 -k2,2n > shH2AZ_total_Input_125bp_control_lambda.sorted.bdg
bedGraphToBigWig shH2AZ_total_Input_125bp_control_lambda.sorted.bdg ~/Data/References/Genomes/Homo_sapiens/GRCh37_hg19_UCSC/chromSizes.txt shH2AZ_total_Input_125bp_control_lambda.sorted.bw
