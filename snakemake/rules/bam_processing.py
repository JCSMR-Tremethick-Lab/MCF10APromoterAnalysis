__author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2016-04-23"

# vim: syntax=python tabstop=4 expandtab
# coding: utf-8


"""
Rules for processing SAM/BAM files
For usage, include this in your workflow.
"""

def bam_merge_input(wildcards):
    fn = []
    for i in config["sample"][wildcards.sample]:
        fn.append("./" + wildcards.processed_dir + "/" + wildcards.genome_version + "/duplicates_removed/" + i + ".DeDup.sorted.fastq_q20.bam")
    return(fn)

wrapper_dir = "/home/skurscheid/Development/snakemake-wrappers/bio"

# import other packages
import os
import fnmatch
from snakemake.exceptions import MissingInputException

# rules
rule all:
    input:
        expand("./processed_data/duplicates_removed/{unit}.DeDup.sorted.fastq_q20.bam", unit = config["units"])

rule bam_merge_dummy:
    input:
        expand("./{processed_dir}/{genome_version}/duplicates_removed/{sample}.Q20.DeDup.sorted.bam", sample = config["sample"], processed_dir = config["processed_dir"], genome_version = "hg38"),
        expand("./{processed_dir}/{genome_version}/duplicates_removed/merged/{sample}.Q20.DeDup.sorted.bam.bai", sample = config["sample"], processed_dir = config["processed_dir"], genome_version = "hg38")


rule bam_sort:
    version:
        "0.1"
    params:
        qual = config["alignment_quality"],
        data_dir = config["data_dir"]
    input:
        "{params.data_dir}/{unit}.fastq_q20.bam"
    output:
        "./processed_data/sorted/{unit}.sorted.fastq_q20.bam"
    shell:
        "samtools sort {input} -T {wildcards.unit}.Q{params.qual}.sorted -o {output}"

rule bam_mark_duplicates:
    version:
        "0.1"
    params:
        qual = config["alignment_quality"]
    input:
        "./processed_data/sorted/{unit}.sorted.fastq_q20.bam"
    output:
        "./processed_data/duplicates_marked/{unit}.MkDup.sorted.fastq_q20.bam"
    shell:
        """
            java -Djava.io.tmpdir=/home/skurscheid/tmp \
            -Xmx36G \
            -jar /home/skurscheid/Bioinformatics/picard-tools-1.131/picard.jar MarkDuplicates \
            INPUT={input}\
            OUTPUT={output}\
            ASSUME_SORTED=TRUE\
            METRICS_FILE={output}.metrics.txt
        """

rule bam_index:
    version:
        "0.1"
    params:
        qual = config["alignment_quality"]
    input:
        "./processed_data/duplicates_marked/{unit}.MkDup.sorted.fastq_q20.bam"
    output:
        "./processed_data/duplicates_marked/{unit}.MkDup.sorted.fastq_q20.bam.bai"
    shell:
        "cd processed_data/duplicates_marked && samtools index ../.{input}"

rule bam_rmdup:
    version:
        "0.1"
    input:
        "./processed_data/duplicates_marked/{unit}.MkDup.sorted.fastq_q20.bam"
    output:
        "./processed_data/duplicates_removed/{unit}.DeDup.sorted.fastq_q20.bam",
        "./processed_data/duplicates_removed/{unit}.DeDup.sorted.fastq_q20.bam.bai"
    shell:
        "samtools rmdup {input} {output[0]}; samtools index {output[0]} {output[1]}"

rule bam_merge:
    version:
        0.1
    input:
        bam_merge_input
    output:
        protected("./{processed_dir}/{genome_version}/duplicates_removed/{sample}.Q20.DeDup.sorted.bam")
    wrapper:
        "file://" + wrapper_dir + "/samtools/merge/wrapper.py"

rule index_merged_bam:
    input:
        rules.bam_merge.output
    output:
        protected("./{processed_dir}/{genome_version}/duplicates_removed/merged/{sample}.Q{qual}.sorted.DeDup.bam.bai")
    shell:
        "samtools index {input} {output}"
