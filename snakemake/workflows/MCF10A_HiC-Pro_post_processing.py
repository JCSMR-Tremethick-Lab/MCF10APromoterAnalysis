_author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2018-02-13"

from snakemake.exceptions import MissingInputException
import os
import glob
from os.path import join

rule:
    version: 0.1

localrules:
    all

home = os.environ['HOME']
wrapper_dir = os.environ['HOME'] + "/Development/snakemake-wrappers/bio"
include_prefix= os.environ['HOME'] + "/Development/JCSMR-Tremethick-Lab/Breast/snakemake/rules/"
armatus_bin = os.environ['HOME'] + "/bin/armatus-linux-x64"
bedtools_location = os.environ['HOME'] + "/miniconda3/envs/chrom3d/bin/"
python_bin = os.environ['HOME'] + "/miniconda3/envs/py27/bin/python"

SAMPLES = ["MCF10A_shH2AZ_1",
           "MCF10A_shH2AZ_2",
           "MCF10AT1_1",
           "MCF10AT1_2",
           "MCF10A_TGFb_1",
           "MCF10A_TGFb_2",
           "MCF10A_WT_1",
           "MCF10A_WT_2",
           "MCF10CA1A_1",
           "MCF10CA1A_2",
           "MCF10CA1A_3"]

rule ice_normalisation:
    version:
        0.1
    params:
        pythonBin = python_bin,
        iceBin = home + "/bin/HiC-Pro_2.10.0/scripts/ice",
        filter_low_counts_perc = "0.02",
        filter_high_counts_perc = "0",
        max_iter = "100",
        eps = "0.1",
        output_bias = "1",
        verbose = "1"
    log:
        "logs/{sample}/ice_{distance}.log"
    input:
        raw_contacts = "hic_results/matrix/{sample}/raw/{distance}/{sample}_{distance}.matrix"
    output:
        iced_contacts = "hic_results/matrix/{sample}/iced/{distance}/{sample}_{distance}_iced.matrix"
    shell:
        """
            {params.pythonBin} {params.iceBin} --results_filename {output.iced_contacts}\
                                               --filter_low_counts_perc {params.filter_low_counts_perc}\
                                               --filter_high_counts_perc {params.filter_high_counts_perc}\
                                               --max_iter {params.max_iter}\
                                               -eps {params.eps}\
                                               --remove-all-zeros-loci\
                                               --output-bias {params.output_bias}\
                                               --verbose {params.verbose}\
                                               {input.raw_contacts} >> {log}
        """

iced_matrices=expand("hic_results/matrix/{sample}/iced/{distance}/{sample}_{distance}_iced.matrix",
                     sample=SAMPLES,
                     distance=["1000000", "150000", "40000", "500000"])


rule convert_list_to_coo:
    version:
        0.1
    params:
        bash_dir = os.environ['HOME'] + "/Development/JCSMR-Tremethick-Lab/Breast/bash"
    input:
        "40k_list/{sample}.{chr1}.{chr2}.{res}.txt"
    output:
        "intra_chr_RAWobserved/{sample}.{chr1}.{chr2}.{res}.coo"
    shell:
    	"""
    		{params.bash_dir}/convert_list_to_coo.sh {input} {output}
    	"""

rule convert_1MB_inter_list_to_coo:
    version:
        0.1
    params:
        bash_dir = os.environ['HOME'] + "/Development/JCSMR-Tremethick-Lab/Breast/bash"
    input:
        "1M_intra_inter_list/{sample}.{chr1}.{chr2}.{res}.txt"
    output:
        "1M_inter_chr_RAWobserved/{sample}.{chr1}.{chr2}.{res}.coo"
    shell:
    	"""
    		{params.bash_dir}/convert_list_to_coo_1M_inter.sh {input} {output}
    	"""

rule make_NCHG_input:
    version:
        0.1
    params:
        bash_dir = os.environ['HOME'] + "/Development/JCSMR-Tremethick-Lab/Breast/bash"
    input:
        coo = "intra_chr_RAWobserved/{sample}.{chr1}.{chr2}.{res}.coo",
        domains = "40k_intra/{sample}.{chr1}.{chr2}.{res}.consensus.domains"
    output:
        "intra_chr_bedpe/{sample}.{chr1}.{chr2}.{res}.domains.RAW.bedpe"
    shell:
        """
            {params.bash_dir}/make_NCHG_input.sh {input.domains} {input.coo} chr{wildcards.chr1} > {output}
        """

rule make_interchr_NCHG_input:
    version:
        0.1
    params:
        python_dir = os.environ['HOME'] + "/Development/JCSMR-Tremethick-Lab/Breast/Python",
        python_bin = os.environ['HOME'] + "/miniconda3/envs/chrom3d/bin/python"
    input:
        coo = "1M_inter_chr_RAWobserved/{sample}.{chr1}.{chr2}.{res}.coo",
        blacklist = os.environ['HOME'] + "/Data/References/Annotations/Homo_sapiens/hg19/UCSC/unmappable_blacklist.bed",
        genomeSizeFile = os.environ['HOME'] + "/Data/References/Annotations/Homo_sapiens/hg19/UCSC/hg19.chrom.sizes.sorted"
    output:
        "inter_chr_bedpe/{sample}.{chr1}.{chr2}.{res}.RAW.bedpe"
    shell:
        """
            {params.python_bin} {params.python_dir}/make_interchr_NCHG_input.py {input.coo} {input.blacklist} {input.genomeSizeFile} chr{wildcards.chr1} chr{wildcards.chr2} > {output}
        """

rule all:
    input:
#        expand("40k_intra/{sample}.domains",
#              sample = SAMPLES),
        # expand("intra_chr_bedpe/{sample}.domains.RAW.bedpe",
        #        sample = SAMPLES),
        # expand("inter_chr_bedpe/{sample}.RAW.bedpe",
        #        sample = INTERCHR)
    	iced_matrices
