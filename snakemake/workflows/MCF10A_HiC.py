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
        "40k_intra/{sample}.{chr1}.{chr2}.{res}.consensus.txt"
    params:
        prefix = "{sample}.{chr1}.{chr2}.{res}"
    shell:
        """
            ~/bin/armatus-linux-x64 -i {input} -g 1.0 -s 0.5 -o {params.prefix} -r 40000 -c {wildcards.chr1}
        """

rule create_domains_file:
    version:
        0.1
    input:
        "40k_intra/{sample}.{chr1}.{chr2}.{res}.consensus.txt"
    output:
        "40k_intra/{sample}.{chr1}.{chr2}.{res}.domains"
    params:
        genomeSizeFile = os.environ['HOME'] + "/Data/References/Annotations/Homo_sapiens/hg19/UCSC/hg19.chrom.sizes.sorted",
        bedtools_location = bedtools_location
    shell:
        """"
            cat {input} | sort -k1,1 -k2,2n - | {params.bedtools_location}/mergeBed -i - | awk '{print "chr"$1 "\\t" $2 "\\t" $3 "\\tdomain"}' - | {params.bedtools_location}/complementBed -i {input} -g {params.genomeSizeFile} | sort -k1,1 -k2,2n - | awk '{if($4=="domain") print $0; else print $1 "\\t" $2 "\\t" $3 "\\tgap"}' > {output}
        """"

rule convert_list_to_coo:
    version:
        0.1
    input:
        "40k_list/{sample}.{chr1}.{chr2}.{res}.txt"
    output:
        "40k_list/{sample}.{chr1}.{chr2}.{res}.coo"
    shell:
        """
            sed 's/-/ /g' {input} | sed 's/k/000/g' | awk '{if($5!=0.0) printf("%s\t%s\t%s\n",$2-40000,$4-40000,$5)}' > {output}
        """

rule all:
    input:
        expand("40k_intra/{sample}.domains",
               sample = SAMPLES),
        expand("40k_list/{sample}.coo",
               sample = SAMPLES)
