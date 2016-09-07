__author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2016-09-07"

# vim: syntax=python tabstop=4 expandtab
# coding: utf-8


"""
Rules for post-processing BAM files
For usage, include this in your workflow.
"""

wrapper_dir = "/home/skurscheid/Development/snakemake-wrappers/bio"

rule target_bedtools_bamtofastq_pe:
    input:
        expand("./fastq/{sample}.{read}.fastq",
                sample = config["units"],
                read = ("end1", "end2"))

rule bedtools_bamtofastq_pe:
    message:
        "Extracting FASTQ reads from BAM file in paired-end mode..."
    params:
    input:
        bam_file = "./bam/{sample}.fastq_q20.bam"
    output:
        read1 = "./fastq/{sample}.end1.fastq",
        read2 = "./fastq/{sample}.end2.fastq"
    wrapper:
        "file://" + wrapper_dir + "/bedtools/bamtofastq/wrapper.py"
