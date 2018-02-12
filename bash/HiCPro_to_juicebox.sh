#!/bin/bash
# for sample MCF10A_WT
cd /home/sebastian/Data/Tremethick/Breast/HiC/processed_data/HiCPro_output/MCF10A_WT

~/Bioinformatics/HiC-Pro_2.9.0/bin/utils/hicpro2juicebox.sh -i hic_results/data/MCF10A_WT_sample1_rep1/MCF10A_WT_sample1_rep1_allValidPairs\
                                                            -g /home/sebastian/Data/References/Genomes/Homo_sapiens/GRCh38_ensembl84/canonical/GRCh38.chromSizes\
                                                            -j /home/sebastian/Bioinformatics/Juicebox/juicer_tools.1.8.9_jcuda.0.8.jar\
                                                            1>>hicpro2juicebox.log 2>>hicpro2juicebox.err &

~/Bioinformatics/HiC-Pro_2.9.0/bin/utils/hicpro2juicebox.sh -i hic_results/data/MCF10A_WT_sample1_rep2/MCF10A_WT_sample1_rep2_allValidPairs\
                                                            -g /home/sebastian/Data/References/Genomes/Homo_sapiens/GRCh38_ensembl84/canonical/GRCh38.chromSizes\
                                                            -j /home/sebastian/Bioinformatics/Juicebox/juicer_tools.1.8.9_jcuda.0.8.jar\
                                                            -t juicebox_data/MCF10A_WT_sample1_rep2/temp \
                                                            -o juicebox_data/MCF10A_WT_sample1_rep2 \
                                                            1>>hicpro2juicebox.log 2>>hicpro2juicebox.err &


~/Bioinformatics/HiC-Pro_2.9.0/bin/utils/hicpro2juicebox.sh -i hic_results/data/MCF10A_WT_sample2_rep1/MCF10A_WT_sample2_rep1_allValidPairs\
                                                            -g /home/sebastian/Data/References/Genomes/Homo_sapiens/GRCh38_ensembl84/canonical/GRCh38.chromSizes\
                                                            -j /home/sebastian/Bioinformatics/Juicebox/juicer_tools.1.8.9_jcuda.0.8.jar\
                                                            -t juicebox_data/MCF10A_WT_sample2_rep1/temp \
                                                            -o juicebox_data/MCF10A_WT_sample2_rep1 \
                                                            1>>hicpro2juicebox.log 2>>hicpro2juicebox.err &


~/Bioinformatics/HiC-Pro_2.9.0/bin/utils/hicpro2juicebox.sh -i hic_results/data/MCF10A_WT_sample2_rep2/MCF10A_WT_sample2_rep2_allValidPairs\
                                                            -g /home/sebastian/Data/References/Genomes/Homo_sapiens/GRCh38_ensembl84/canonical/GRCh38.chromSizes\
                                                            -j /home/sebastian/Bioinformatics/Juicebox/juicer_tools.1.8.9_jcuda.0.8.jar\
                                                            -t juicebox_data/MCF10A_WT_sample2_rep2/temp \
                                                            -o juicebox_data/MCF10A_WT_sample2_rep2 \
                                                            1>>hicpro2juicebox.log 2>>hicpro2juicebox.err &

######################################################################################
cd /home/sebastian/Data/Tremethick/Breast/HiC/processed_data/HiCPro_output/MCF10A_TGFb

export sample="MCF10A_TGFb_sample1_rep1"
~/Bioinformatics/HiC-Pro_2.9.0/bin/utils/hicpro2juicebox.sh -i hic_results/data/${sample}/${sample}_allValidPairs\
                                                            -g /home/sebastian/Data/References/Genomes/Homo_sapiens/GRCh38_ensembl84/canonical/GRCh38.chromSizes\
                                                            -j /home/sebastian/Bioinformatics/Juicebox/juicer_tools.1.8.9_jcuda.0.8.jar\
                                                            -t juicebox_data/${sample}/temp \
                                                            -o juicebox_data/${sample} \
                                                            1>>hicpro2juicebox.log 2>>hicpro2juicebox.err &

export sample="MCF10A_TGFb_sample1_rep2"
~/Bioinformatics/HiC-Pro_2.9.0/bin/utils/hicpro2juicebox.sh -i hic_results/data/${sample}/${sample}_allValidPairs\
                                                            -g /home/sebastian/Data/References/Genomes/Homo_sapiens/GRCh38_ensembl84/canonical/GRCh38.chromSizes\
                                                            -j /home/sebastian/Bioinformatics/Juicebox/juicer_tools.1.8.9_jcuda.0.8.jar\
                                                            -t juicebox_data/${sample}/temp \
                                                            -o juicebox_data/${sample} \
                                                            1>>hicpro2juicebox.log 2>>hicpro2juicebox.err &


export sample="MCF10A_TGFb_sample2_rep1"
~/Bioinformatics/HiC-Pro_2.9.0/bin/utils/hicpro2juicebox.sh -i hic_results/data/${sample}/${sample}_allValidPairs\
                                                            -g /home/sebastian/Data/References/Genomes/Homo_sapiens/GRCh38_ensembl84/canonical/GRCh38.chromSizes\
                                                            -j /home/sebastian/Bioinformatics/Juicebox/juicer_tools.1.8.9_jcuda.0.8.jar\
                                                            -t juicebox_data/${sample}/temp \
                                                            -o juicebox_data/${sample} \
                                                            1>>hicpro2juicebox.log 2>>hicpro2juicebox.err &


export sample="MCF10A_TGFb_sample2_rep2"
~/Bioinformatics/HiC-Pro_2.9.0/bin/utils/hicpro2juicebox.sh -i hic_results/data/${sample}/${sample}_allValidPairs\
                                                            -g /home/sebastian/Data/References/Genomes/Homo_sapiens/GRCh38_ensembl84/canonical/GRCh38.chromSizes\
                                                            -j /home/sebastian/Bioinformatics/Juicebox/juicer_tools.1.8.9_jcuda.0.8.jar\
                                                            -t juicebox_data/${sample}/temp \
                                                            -o juicebox_data/${sample} \
                                                            1>>hicpro2juicebox.log 2>>hicpro2juicebox.err &

######################################################################################
cd /home/sebastian/Data/Tremethick/Breast/HiC/processed_data/HiCPro_output/MCF10CA1A

export sample="MCF10CA1A_sample1_rep1"
~/Bioinformatics/HiC-Pro_2.9.0/bin/utils/hicpro2juicebox.sh -i hic_results/data/${sample}/${sample}_allValidPairs\
                                                            -g /home/sebastian/Data/References/Genomes/Homo_sapiens/GRCh38_ensembl84/canonical/GRCh38.chromSizes\
                                                            -j /home/sebastian/Bioinformatics/Juicebox/juicer_tools.1.8.9_jcuda.0.8.jar\
                                                            -t juicebox_data/${sample}/temp \
                                                            -o juicebox_data/${sample} \
                                                            1>>hicpro2juicebox.log 2>>hicpro2juicebox.err &

export sample="MCF10CA1A_sample1_rep2"
~/Bioinformatics/HiC-Pro_2.9.0/bin/utils/hicpro2juicebox.sh -i hic_results/data/${sample}/${sample}_allValidPairs\
                                                            -g /home/sebastian/Data/References/Genomes/Homo_sapiens/GRCh38_ensembl84/canonical/GRCh38.chromSizes\
                                                            -j /home/sebastian/Bioinformatics/Juicebox/juicer_tools.1.8.9_jcuda.0.8.jar\
                                                            -t juicebox_data/${sample}/temp \
                                                            -o juicebox_data/${sample} \
                                                            1>>hicpro2juicebox.log 2>>hicpro2juicebox.err &


export sample="MCF10CA1A_sample2_rep1"
~/Bioinformatics/HiC-Pro_2.9.0/bin/utils/hicpro2juicebox.sh -i hic_results/data/${sample}/${sample}_allValidPairs\
                                                            -g /home/sebastian/Data/References/Genomes/Homo_sapiens/GRCh38_ensembl84/canonical/GRCh38.chromSizes\
                                                            -j /home/sebastian/Bioinformatics/Juicebox/juicer_tools.1.8.9_jcuda.0.8.jar\
                                                            -t juicebox_data/${sample}/temp \
                                                            -o juicebox_data/${sample} \
                                                            1>>hicpro2juicebox.log 2>>hicpro2juicebox.err &


export sample="MCF10CA1A_sample2_rep2"
~/Bioinformatics/HiC-Pro_2.9.0/bin/utils/hicpro2juicebox.sh -i hic_results/data/${sample}/${sample}_allValidPairs\
                                                            -g /home/sebastian/Data/References/Genomes/Homo_sapiens/GRCh38_ensembl84/canonical/GRCh38.chromSizes\
                                                            -j /home/sebastian/Bioinformatics/Juicebox/juicer_tools.1.8.9_jcuda.0.8.jar\
                                                            -t juicebox_data/${sample}/temp \
                                                            -o juicebox_data/${sample} \
                                                            1>>hicpro2juicebox.log 2>>hicpro2juicebox.err &

export sample="MCF10CA1A_sample3_rep1"
~/Bioinformatics/HiC-Pro_2.9.0/bin/utils/hicpro2juicebox.sh -i hic_results/data/${sample}/${sample}_allValidPairs\
                                                            -g /home/sebastian/Data/References/Genomes/Homo_sapiens/GRCh38_ensembl84/canonical/GRCh38.chromSizes\
                                                            -j /home/sebastian/Bioinformatics/Juicebox/juicer_tools.1.8.9_jcuda.0.8.jar\
                                                            -t juicebox_data/${sample}/temp \
                                                            -o juicebox_data/${sample} \
                                                            1>>hicpro2juicebox.log 2>>hicpro2juicebox.err &


export sample="MCF10CA1A_sample3_rep2"
~/Bioinformatics/HiC-Pro_2.9.0/bin/utils/hicpro2juicebox.sh -i hic_results/data/${sample}/${sample}_allValidPairs\
                                                            -g /home/sebastian/Data/References/Genomes/Homo_sapiens/GRCh38_ensembl84/canonical/GRCh38.chromSizes\
                                                            -j /home/sebastian/Bioinformatics/Juicebox/juicer_tools.1.8.9_jcuda.0.8.jar\
                                                            -t juicebox_data/${sample}/temp \
                                                            -o juicebox_data/${sample} \
                                                            1>>hicpro2juicebox.log 2>>hicpro2juicebox.err &
