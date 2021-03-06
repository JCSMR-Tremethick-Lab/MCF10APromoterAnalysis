__author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2017-02-01"

# vim: syntax=python tabstop=4 expandtab
# coding: utf-8


"""
Rules for processing SAM/BAM files
For usage, include this in your workflow.
"""

# import other packages
import os
import fnmatch
from snakemake.exceptions import MissingInputException

# set some local variables
home = os.environ['HOME']
REF_GENOME = config["references"]["genomes"][1]

def bam_merge_input(wildcards):
    fn = []
    path = "/".join((wildcards["assayID"],
                     wildcards["runID"],
                     wildcards["outdir"],
                     wildcards["reference_version"],
                     "bowtie2",
                     wildcards["duplicates"]))
    for i in config["samples"]["ChIP-Seq"]["replicates"][wildcards["sampleGroup"]]:
        fn.append("/".join((path, ".".join((i, ("".join("Q", config["alignment_quality"])),"sorted.bam")))))
    return(fn)

rule run_bam_merge:
    input:
        expand("{assayID}/{runID}/{outdir}/{reference_version}/{application}/{command}/{duplicates}/{sampleGroup}.{suffix}",
               assayID = "ChIP-Seq",
               runID = "merged",
               outdir = config["processed_dir"],
               reference_version = config["references"][REF_GENOME]["version"][0],
               application = "samtools",
               command = "merge",
               duplicates = ["duplicates_marked", "duplicates_removed"],
               sampleGroup = ["H2AZ_10A_high", "H2AZ_TGFb_10A", "Inp_10A_WT_high", "Inp_10A_TGFb_high", "Inp_shZ_10A_high"],
               suffix = ["bam", "bam.bai"])

rule bam_merge:
    version:
        0.3
    params:
        cwd = os.getcwd()
    threads:
        lambda wildcards: int(str(config["program_parameters"]["bt2_params"]["threads"]).strip("['']"))
    input:
        lambda wildcards: expand("{assayID}/{runID}/{outdir}/{reference_version}/bowtie2/{duplicates}/{unit}.Q{qual}.{suffix}",
                                 assayID = wildcards["assayID"],
                                 runID = wildcards["runID"],
                                 outdir = wildcards["outdir"],
                                 reference_version = wildcards["reference_version"],
                                 duplicates = wildcards["duplicates"],
                                 unit = config["samples"]["ChIP-Seq"]["replicates"][wildcards["sampleGroup"]],
                                 qual = config["alignment_quality"],
                                 suffix = "sorted.bam")
    output:
        protected("{assayID}/{runID}/{outdir}/{reference_version}/{application}/{command}/{duplicates}/{sampleGroup}.bam")
    run:
        if (len(input) > 1):
            shell("samtools merge --threads {threads} {output} {input}")
        else:
            shell("ln -s {params.cwd}/{input} {output}")

rule index_merged_bam:
    input:
        rules.bam_merge.output
    output:
        protected("{assayID}/{runID}/{outdir}/{reference_version}/{application}/{command}/{duplicates}/{sampleGroup}.bam.bai")
    shell:
        "samtools index {input} {output}"
