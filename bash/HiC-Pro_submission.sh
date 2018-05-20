# step wise processing with intermediate saving
export RES_PREFIX="/home/skurscheid/Data/Tremethick/Breast/HiC/HiCPro_output_run3"
export INPUT="/home/skurscheid/Data/Tremethick/Breast/HiC/HiCPro_input_run2"

~/bin/HiC-Pro_2.10.0/bin/HiC-Pro -i ${INPUT} -o ${RES_PREFIX} -c ${INPUT}/config-hicpro.txt -s mapping -s quality_checks -p
cd ${RES_PREFIX}; qsub HiCPro_step1_hicpro.sh
qsub HiCPro_step2_hicpro.sh ; cd ..
~/bin/HiC-Pro_2.10.0/bin/HiC-Pro -i ${RES_PREFIX}/bowtie_results/bwt2 -o ${RES_PREFIX}.1/ -c ${INPUT}/config-hicpro.txt -s proc_hic -s quality_checks -p
cd ${RES_PREFIX}.1; qsub HiCPro_step1_hicpro.sh
qsub HiCPro_step2_hicpro.sh ; cd ..
~/bin/HiC-Pro_2.10.0/bin/HiC-Pro -i ${RES_PREFIX}.1/hic_results/data -o ${RES_PREFIX}.2/ -c ${INPUT}/config-hicpro.txt -s build_contact_maps -s ice_norm -p
cd ${RES_PREFIX}.2; qsub HiCPro_step1_hicpro.sh
qsub HiCPro_step2_hicpro.sh ; cd ..
