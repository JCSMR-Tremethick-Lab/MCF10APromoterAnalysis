__author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2016-09-07"

# vim: syntax=python tabstop=4 expandtab
# coding: utf-8


"""
Rules for running MACS2 for ChIP peak calling.
For usage, include this in your workflow.
"""

def get_replicates_input(wildcards):
    fn = []
    for i in config["samples"][wildcards.digest][wildcards.Input][wildcards.sample + "_" + wildcards.Input + "_" + wildcards.digest]:
        fn.append("processed_data/hg38/duplicates_removed/" + i + ".DeDup.sorted.fastq_q20.bam")
    return(fn)

def get_replicates_chip(wildcards):
    fn = []
    for i in config["samples"][wildcards.digest][wildcards.ChIP][wildcards.sample + "_" + wildcards.ChIP + "_" + wildcards.digest]:
        fn.append("processed_data/hg38/duplicates_removed/" + i + ".DeDup.sorted.fastq_q20.bam")
    return(fn)


rule macs2_callpeak_dummy:
    input:
        expand("processed_data/hg38/macs2/callpeak/{digest}/{ChIP}_vs_{Input}/{sample}",
               digest = "H",
               ChIP = "H2AZ",
               Input = "Input",
               sample = config["samples"]["sample"])

rule macs2_callpeak:
    params:
        extsize = config["parameters"]["macs2"]["extsize"],
        macs2_dir = config["macs2_dir"]
    input:
        input = get_replicates_input,
        chip = get_replicates_chip
    output:
        "processed_data/hg38/macs2/callpeak/{digest}/{ChIP}_vs_{Input}/{sample}"
    shell:
        """
            {params.macs2_dir}/macs2 callpeak -B \
                                              -t {input.chip}\
                                              -c {input.input}\
                                              -n {wildcards.sample}\
                                              --nomodel\
                                              --extsize {params.extsize}\
                                              --outdir {output}
        """
