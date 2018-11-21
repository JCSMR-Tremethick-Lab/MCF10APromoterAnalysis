export dataDir="/home/sebastian/Data/Tremethick/Breast/PromoterSeqCap/processed_data/GRCh37_hg19_UCSC/bowtie2"
export outputDir="/home/sebastian/Data/Tremethick/Breast/PromoterSeqCap/processed_data/GRCh37_hg19_UCSC/macs2/callpeak_nocutoff"

~/miniconda3/envs/py27/bin/macs2 callpeak -g 40999507\
                                          -t ${dataDir}/A_H2AZ_H_r1_R1.final.bam\
                                          -c ${dataDir}/A_Inp_H_r1_R1.final.bam\
                                             ${dataDir}/A_Inp_H_r2_R1.final.bam\
                                          --qvalue 0.99\
                                          --outdir ${outputDir}/A_H2AZ_H_r1_R1\
                                          -n A_H2AZ_H_r1_R1\
                                          --bdg\
                                          --SPMR\
                                          --cutoff-analysis\
                                          --call-summits\
                                          --nomodel\
                                          --extsize 147 &

~/miniconda3/envs/py27/bin/macs2 callpeak -g 40999507\
                                          -t ${dataDir}/A_H2AZ_H_r2_R1.final.bam\
                                          -c ${dataDir}/A_Inp_H_r1_R1.final.bam\
                                             ${dataDir}/A_Inp_H_r2_R1.final.bam\
                                          --qvalue 0.99\
                                          --outdir ${outputDir}/A_H2AZ_H_r2_R1\
                                          -n A_H2AZ_H_r2_R1\
                                          --bdg\
                                          --SPMR\
                                          --cutoff-analysis\
                                          --call-summits\
                                          --nomodel\
                                          --extsize 147 &

~/miniconda3/envs/py27/bin/macs2 callpeak -g 40999507\
                                          -t ${dataDir}/A_H2AZ_L_r1_R1.final.bam\
                                          -c ${dataDir}/A_Inp_L_r1_R1.final.bam\
                                             ${dataDir}/A_Inp_L_r2_R1.final.bam\
                                          --outdir ${outputDir}/A_H2AZ_L_r1_R1\
                                          -n A_H2AZ_L_r1_R1\
                                          --trackline\
                                          --SPMR\
                                          --cutoff-analysis\
                                          --call-summits\
                                          --nomodel\
                                          --extsize 147 &

~/miniconda3/envs/py27/bin/macs2 callpeak -g 40999507\
                                          -t ${dataDir}/A_H2AZ_L_r2_R1.final.bam\
                                          -c ${dataDir}/A_Inp_L_r1_R1.final.bam\
                                             ${dataDir}/A_Inp_L_r2_R1.final.bam\
                                          --outdir ${outputDir}/A_H2AZ_L_r2_R1 \
                                          -n A_H2AZ_L_r2_R1\
                                          --trackline\
                                          --SPMR\
                                          --cutoff-analysis\
                                          --call-summits\
                                          --nomodel\
                                          --extsize 147 &

# running IDR
~/miniconda3/envs/chip-seq/bin/idr --samples ${outputDir}/A_H2AZ_H_r1_R1/A_H2AZ_H_r1_R1_peaks.narrowPeak ${outputDir}/A_H2AZ_H_r2_R1/A_H2AZ_H_r2_R1_peaks.narrowPeak \
                                   --input-file-type narrowPeak\
                                   --output-file /home/sebastian/Data/Tremethick/Breast/PromoterSeqCap/processed_data/GRCh37_hg19_UCSC/idr_noMacs2Cutoff/H2AZ_H_idr_values.txt \
                                   --log-output-file /home/sebastian/Data/Tremethick/Breast/PromoterSeqCap/processed_data/GRCh37_hg19_UCSC/idr_noMacs2Cutoff/H2AZ_H_idr_log.txt \
                                   --soft-idr-threshold 0.1 \
                                   --plot
