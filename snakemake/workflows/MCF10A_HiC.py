_author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2017-08-02"

from snakemake.exceptions import MissingInputException
import os
import glob

rule:
    version: 0.1

localrules:
    all

wrapper_dir = os.environ['HOME'] + "/Development/snakemake-wrappers/bio"
include_prefix= os.environ['HOME'] + "/Development/JCSMR-Tremethick-Lab/Breast/snakemake/rules/"
armatus_bin = os.environ['HOME'] + "/bin/armatus-linux-x64"
bedtools_location = os.environ['HOME'] + "/miniconda3/envs/chrom3d/bin/"

SAMPLES = ["shZ-rep2.10.10.40k" ,
           "shZ-rep2.11.11.40k" ,
           "shZ-rep2.1.1.40k" ,
           "shZ-rep2.12.12.40k" ,
           "shZ-rep2.13.13.40k" ,
           "shZ-rep2.14.14.40k" ,
           "shZ-rep2.15.15.40k" ,
           "shZ-rep2.16.16.40k" ,
           "shZ-rep2.17.17.40k" ,
           "shZ-rep2.18.18.40k" ,
           "shZ-rep2.19.19.40k" ,
           "shZ-rep2.20.20.40k" ,
           "shZ-rep2.21.21.40k" ,
           "shZ-rep2.22.22.40k" ,
           "shZ-rep2.2.2.40k" ,
           "shZ-rep2.3.3.40k" ,
           "shZ-rep2.4.4.40k" ,
           "shZ-rep2.5.5.40k" ,
           "shZ-rep2.6.6.40k" ,
           "shZ-rep2.7.7.40k" ,
           "shZ-rep2.8.8.40k" ,
           "shZ-rep2.9.9.40k" ,
           "shZ-rep2.X.X.40k"]

INTERCHR = ["shZ-rep2.10.11.1000k",
            "shZ-rep2.10.12.1000k",
            "shZ-rep2.10.13.1000k",
            "shZ-rep2.10.14.1000k",
            "shZ-rep2.10.15.1000k",
            "shZ-rep2.10.16.1000k",
            "shZ-rep2.10.17.1000k",
            "shZ-rep2.10.18.1000k",
            "shZ-rep2.10.19.1000k",
            "shZ-rep2.10.20.1000k",
            "shZ-rep2.10.21.1000k",
            "shZ-rep2.10.22.1000k",
            "shZ-rep2.10.X.1000k",
            "shZ-rep2.1.10.1000k",
            "shZ-rep2.1.11.1000k",
            "shZ-rep2.11.12.1000k",
            "shZ-rep2.11.13.1000k",
            "shZ-rep2.11.14.1000k",
            "shZ-rep2.11.15.1000k",
            "shZ-rep2.11.16.1000k",
            "shZ-rep2.11.17.1000k",
            "shZ-rep2.11.18.1000k",
            "shZ-rep2.11.19.1000k",
            "shZ-rep2.11.20.1000k",
            "shZ-rep2.1.12.1000k",
            "shZ-rep2.11.21.1000k",
            "shZ-rep2.11.22.1000k",
            "shZ-rep2.1.13.1000k",
            "shZ-rep2.1.14.1000k",
            "shZ-rep2.1.15.1000k",
            "shZ-rep2.1.16.1000k",
            "shZ-rep2.1.17.1000k",
            "shZ-rep2.1.18.1000k",
            "shZ-rep2.1.19.1000k",
            "shZ-rep2.11.X.1000k",
            "shZ-rep2.1.20.1000k",
            "shZ-rep2.1.2.1000k",
            "shZ-rep2.1.21.1000k",
            "shZ-rep2.12.13.1000k",
            "shZ-rep2.12.14.1000k",
            "shZ-rep2.12.15.1000k",
            "shZ-rep2.12.16.1000k",
            "shZ-rep2.12.17.1000k",
            "shZ-rep2.12.18.1000k",
            "shZ-rep2.12.19.1000k",
            "shZ-rep2.12.20.1000k",
            "shZ-rep2.1.22.1000k",
            "shZ-rep2.12.21.1000k",
            "shZ-rep2.12.22.1000k",
            "shZ-rep2.12.X.1000k",
            "shZ-rep2.1.3.1000k",
            "shZ-rep2.13.14.1000k",
            "shZ-rep2.13.15.1000k",
            "shZ-rep2.13.16.1000k",
            "shZ-rep2.13.17.1000k",
            "shZ-rep2.13.18.1000k",
            "shZ-rep2.13.19.1000k",
            "shZ-rep2.13.20.1000k",
            "shZ-rep2.13.21.1000k",
            "shZ-rep2.13.22.1000k",
            "shZ-rep2.13.X.1000k",
            "shZ-rep2.1.4.1000k",
            "shZ-rep2.14.15.1000k",
            "shZ-rep2.14.16.1000k",
            "shZ-rep2.14.17.1000k",
            "shZ-rep2.14.18.1000k",
            "shZ-rep2.14.19.1000k",
            "shZ-rep2.14.20.1000k",
            "shZ-rep2.14.21.1000k",
            "shZ-rep2.14.22.1000k",
            "shZ-rep2.14.X.1000k",
            "shZ-rep2.1.5.1000k",
            "shZ-rep2.15.16.1000k",
            "shZ-rep2.15.17.1000k",
            "shZ-rep2.15.18.1000k",
            "shZ-rep2.15.19.1000k",
            "shZ-rep2.15.20.1000k",
            "shZ-rep2.15.21.1000k",
            "shZ-rep2.15.22.1000k",
            "shZ-rep2.15.X.1000k",
            "shZ-rep2.1.6.1000k",
            "shZ-rep2.16.17.1000k",
            "shZ-rep2.16.18.1000k",
            "shZ-rep2.16.19.1000k",
            "shZ-rep2.16.20.1000k",
            "shZ-rep2.16.21.1000k",
            "shZ-rep2.16.22.1000k",
            "shZ-rep2.16.X.1000k",
            "shZ-rep2.1.7.1000k",
            "shZ-rep2.17.18.1000k",
            "shZ-rep2.17.19.1000k",
            "shZ-rep2.17.20.1000k",
            "shZ-rep2.17.21.1000k",
            "shZ-rep2.17.22.1000k",
            "shZ-rep2.17.X.1000k",
            "shZ-rep2.1.8.1000k",
            "shZ-rep2.18.19.1000k",
            "shZ-rep2.18.20.1000k",
            "shZ-rep2.18.21.1000k",
            "shZ-rep2.18.22.1000k",
            "shZ-rep2.18.X.1000k",
            "shZ-rep2.1.9.1000k",
            "shZ-rep2.19.20.1000k",
            "shZ-rep2.19.21.1000k",
            "shZ-rep2.19.22.1000k",
            "shZ-rep2.19.X.1000k",
            "shZ-rep2.1.X.1000k",
            "shZ-rep2.20.21.1000k",
            "shZ-rep2.20.22.1000k",
            "shZ-rep2.20.X.1000k",
            "shZ-rep2.2.10.1000k",
            "shZ-rep2.2.11.1000k",
            "shZ-rep2.2.12.1000k",
            "shZ-rep2.21.22.1000k",
            "shZ-rep2.2.13.1000k",
            "shZ-rep2.2.14.1000k",
            "shZ-rep2.2.15.1000k",
            "shZ-rep2.2.16.1000k",
            "shZ-rep2.2.17.1000k",
            "shZ-rep2.2.18.1000k",
            "shZ-rep2.2.19.1000k",
            "shZ-rep2.21.X.1000k",
            "shZ-rep2.2.20.1000k",
            "shZ-rep2.2.21.1000k",
            "shZ-rep2.2.22.1000k",
            "shZ-rep2.22.X.1000k",
            "shZ-rep2.2.3.1000k",
            "shZ-rep2.2.4.1000k",
            "shZ-rep2.2.5.1000k",
            "shZ-rep2.2.6.1000k",
            "shZ-rep2.2.7.1000k",
            "shZ-rep2.2.8.1000k",
            "shZ-rep2.2.9.1000k",
            "shZ-rep2.2.X.1000k",
            "shZ-rep2.3.10.1000k",
            "shZ-rep2.3.11.1000k",
            "shZ-rep2.3.12.1000k",
            "shZ-rep2.3.13.1000k",
            "shZ-rep2.3.14.1000k",
            "shZ-rep2.3.15.1000k",
            "shZ-rep2.3.16.1000k",
            "shZ-rep2.3.17.1000k",
            "shZ-rep2.3.18.1000k",
            "shZ-rep2.3.19.1000k",
            "shZ-rep2.3.20.1000k",
            "shZ-rep2.3.21.1000k",
            "shZ-rep2.3.22.1000k",
            "shZ-rep2.3.4.1000k",
            "shZ-rep2.3.5.1000k",
            "shZ-rep2.3.6.1000k",
            "shZ-rep2.3.7.1000k",
            "shZ-rep2.3.8.1000k",
            "shZ-rep2.3.9.1000k",
            "shZ-rep2.3.X.1000k",
            "shZ-rep2.4.10.1000k",
            "shZ-rep2.4.11.1000k",
            "shZ-rep2.4.12.1000k",
            "shZ-rep2.4.13.1000k",
            "shZ-rep2.4.14.1000k",
            "shZ-rep2.4.15.1000k",
            "shZ-rep2.4.16.1000k",
            "shZ-rep2.4.17.1000k",
            "shZ-rep2.4.18.1000k",
            "shZ-rep2.4.19.1000k",
            "shZ-rep2.4.20.1000k",
            "shZ-rep2.4.21.1000k",
            "shZ-rep2.4.22.1000k",
            "shZ-rep2.4.5.1000k",
            "shZ-rep2.4.6.1000k",
            "shZ-rep2.4.7.1000k",
            "shZ-rep2.4.8.1000k",
            "shZ-rep2.4.9.1000k",
            "shZ-rep2.4.X.1000k",
            "shZ-rep2.5.10.1000k",
            "shZ-rep2.5.11.1000k",
            "shZ-rep2.5.12.1000k",
            "shZ-rep2.5.13.1000k",
            "shZ-rep2.5.14.1000k",
            "shZ-rep2.5.15.1000k",
            "shZ-rep2.5.16.1000k",
            "shZ-rep2.5.17.1000k",
            "shZ-rep2.5.18.1000k",
            "shZ-rep2.5.19.1000k",
            "shZ-rep2.5.20.1000k",
            "shZ-rep2.5.21.1000k",
            "shZ-rep2.5.22.1000k",
            "shZ-rep2.5.6.1000k",
            "shZ-rep2.5.7.1000k",
            "shZ-rep2.5.8.1000k",
            "shZ-rep2.5.9.1000k",
            "shZ-rep2.5.X.1000k",
            "shZ-rep2.6.10.1000k",
            "shZ-rep2.6.11.1000k",
            "shZ-rep2.6.12.1000k",
            "shZ-rep2.6.13.1000k",
            "shZ-rep2.6.14.1000k",
            "shZ-rep2.6.15.1000k",
            "shZ-rep2.6.16.1000k",
            "shZ-rep2.6.17.1000k",
            "shZ-rep2.6.18.1000k",
            "shZ-rep2.6.19.1000k",
            "shZ-rep2.6.20.1000k",
            "shZ-rep2.6.21.1000k",
            "shZ-rep2.6.22.1000k",
            "shZ-rep2.6.7.1000k",
            "shZ-rep2.6.8.1000k",
            "shZ-rep2.6.9.1000k",
            "shZ-rep2.6.X.1000k",
            "shZ-rep2.7.10.1000k",
            "shZ-rep2.7.11.1000k",
            "shZ-rep2.7.12.1000k",
            "shZ-rep2.7.13.1000k",
            "shZ-rep2.7.14.1000k",
            "shZ-rep2.7.15.1000k",
            "shZ-rep2.7.16.1000k",
            "shZ-rep2.7.17.1000k",
            "shZ-rep2.7.18.1000k",
            "shZ-rep2.7.19.1000k",
            "shZ-rep2.7.20.1000k",
            "shZ-rep2.7.21.1000k",
            "shZ-rep2.7.22.1000k",
            "shZ-rep2.7.8.1000k",
            "shZ-rep2.7.9.1000k",
            "shZ-rep2.7.X.1000k",
            "shZ-rep2.8.10.1000k",
            "shZ-rep2.8.11.1000k",
            "shZ-rep2.8.12.1000k",
            "shZ-rep2.8.13.1000k",
            "shZ-rep2.8.14.1000k",
            "shZ-rep2.8.15.1000k",
            "shZ-rep2.8.16.1000k",
            "shZ-rep2.8.17.1000k",
            "shZ-rep2.8.18.1000k",
            "shZ-rep2.8.19.1000k",
            "shZ-rep2.8.20.1000k",
            "shZ-rep2.8.21.1000k",
            "shZ-rep2.8.22.1000k",
            "shZ-rep2.8.9.1000k",
            "shZ-rep2.8.X.1000k",
            "shZ-rep2.9.10.1000k",
            "shZ-rep2.9.11.1000k",
            "shZ-rep2.9.12.1000k",
            "shZ-rep2.9.13.1000k",
            "shZ-rep2.9.14.1000k",
            "shZ-rep2.9.15.1000k",
            "shZ-rep2.9.16.1000k",
            "shZ-rep2.9.17.1000k",
            "shZ-rep2.9.18.1000k",
            "shZ-rep2.9.19.1000k",
            "shZ-rep2.9.20.1000k",
            "shZ-rep2.9.21.1000k",
            "shZ-rep2.9.22.1000k",
            "shZ-rep2.9.X.1000k"]

rule prepare_armatus_input:
    version: 0.1
    input:
        "40k_intra/{sample}.{chr1}.{chr2}.{res}.mat"
    output:
        "40k_intra/{sample}.{chr1}.{chr2}.{res}.gz"
    threads: 1
    shell:
        """
            cut -f 2- -d " " {input} | sed "1d" - | tr " " "\t" | gzip > {output}
        """

rule run_armatus:
    version:
        0.1
    input:
        "40k_intra/{sample}.{chr1}.{chr2}.{res}.gz"
    output:
        "40k_intra/{sample}.{chr1}.{chr2}.{res}.consensus"
    params:
        prefix = "{sample}.{chr1}.{chr2}.{res}"
    shell:
        """
            ~/bin/armatus-linux-x64 -i {input} -g 1.0 -s 0.5 -o {params.prefix} -r 40000 -c {wildcards.chr1}
        """

# rule create_domains_file:
#     version:
#         0.1
#     input:
#         "40k_intra/{sample}.{chr1}.{chr2}.{res}.consensus"
#     output:
#         "40k_intra/{sample}.{chr1}.{chr2}.{res}.domains"
#     params:
#         genomeSizeFile = os.environ['HOME'] + "/Data/References/Annotations/Homo_sapiens/hg19/UCSC/hg19.chrom.sizes.sorted",
#         bedtools_location = bedtools_location
#     shell:
#         """"
#         """"

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
        blacklist = "~/Data/References/Annotations/Homo_sapiens/hg19/UCSC/unmappable_blacklist.bed",
        genomeSizeFile = "~/Data/References/Annotations/Homo_sapiens/hg19/UCSC/hg19.chrom.sizes.sorted"
    output:
        "inter_chr_bedpe/{sample}.{chr1}.{chr2}.{res}.domains.RAW.bedpe"
    shell:
        """
            {params.python_bin} {params.python_dir}/make_interchr_NCHG_input.sh {input.coo} {input.blacklist} {input.genomeSizeFile} chr{wildcards.chr1} chr{wildcards.chr2} > {output}
        """

rule all:
    input:
#        expand("40k_intra/{sample}.domains",
#              sample = SAMPLES),
        expand("intra_chr_bedpe/{sample}.domains.RAW.bedpe",
               sample = SAMPLES),
        expand("1M_inter_chr_RAWobserved/{sample}.coo",
               sample = INTERCHR)
