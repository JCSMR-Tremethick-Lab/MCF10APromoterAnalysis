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

prepare_armatus_input:
    input:
        glob_wildcards("{sample}.{chr1}.{chr2}.{res}.mat")
    output:
        "{sample}.{chr1}.{chr2}.{res}.tsv"
    params:
    threads: 1
    shell:
        """
            cut -f 2- -d " " {input} | sed "1d" - | tr " " "\t" > {output}
        """

all:
    input:
        "shZ-rep2.10.10.40k.tsv"
