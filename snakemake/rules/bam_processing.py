__author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2016-02-29"

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

def bam_merge_input(wildcards):
    fn = []
    path = "/".join((wildcards["assayID"],
                     wildcards["runID"],
                     wildcards["outdir"],
                     wildcards["reference_version"],
                     wildcards["application"],
                     wildcards["duplicates"]))
    for i in config["samples"]["ChIP-Seq"]["replicates"][wildcards["sample"]]:
        fn.append("/".join((path, ".".join((i, "Q20.sorted.bam")))))
    return(fn)

rule:
    version: 0.2

# rules
rule bam_quality_filter:
    params:
        qual = config["alignment_quality"]
    input:
        rules.bowtie2_pe.output
    output:
        temp("{assayID}/{runID}/{outdir}/{reference_version}/bowtie2/quality_filtered/{unit}.Q{qual}.bam")
    shell:
        "samtools view -b -h -q {params.qual} {input} > {output}"

rule bam_sort:
    params:
        qual = config["alignment_quality"],
        threads = "2"
    input:
        rules.bam_quality_filter.output
    output:
        "{assayID}/{runID}/{outdir}/{reference_version}/bowtie2/sorted/{unit}.Q{qual}.sorted.bam"
    shell:
        "samtools sort -@ {params.threads} {input} -T {wildcards.unit}.Q{params.qual}.sorted -o {output}"

rule bam_mark_duplicates:
    params:
        qual = config["alignment_quality"],
        picard = home + config["picard"],
        temp = home + config["temp_dir"]
    input:
        rules.bam_sort.output
    output:
        protected("{assayID}/{runID}/{outdir}/{reference_version}/bowtie2/duplicates_marked/{unit}.Q{qual}.sorted.MkDup.bam")
    shell:
        """
            java -Djava.io.tmpdir={params.temp} \
            -Xmx24G \
            -jar {params.picard} MarkDuplicates \
            INPUT={input}\
            OUTPUT={output}\
            ASSUME_SORTED=TRUE\
            METRICS_FILE={output}.metrics.txt
        """

rule bam_index:
    params:
        qual = config["alignment_quality"]
    input:
        rules.bam_mark_duplicates.output
    output:
        protected("{assayID}/{runID}/{outdir}/{reference_version}/bowtie2/duplicates_marked/{unit}.Q{qual}.sorted.MkDup.bam.bai")
    shell:
        "samtools index {input} {output}"

rule bam_rmdup:
    input:
        rules.bam_mark_duplicates.output
    output:
        protected("{assayID}/{runID}/{outdir}/{reference_version}/bowtie2/duplicates_removed/{unit}.Q{qual}.sorted.DeDup.bam")
    shell:
        "samtools rmdup {input} {output}"

rule bam_rmdup_index:
    params:
        qual = config["alignment_quality"]
    input:
        rules.bam_rmdup.output
    output:
        protected("{assayID}/{runID}/{outdir}/{reference_version}/bowtie2/duplicates_removed/{unit}.Q{qual}.sorted.DeDup.bam.bai")
    shell:
        "samtools index {input} {output}"

rule bam_merge:
    version:
        0.2
    input:
        bam_merge_input
    output:
        protected("{assayID}/{runID}/{outdir}/{reference_version}/{application}/{duplicates}/{sample}.bam")
    run:
        if len(input > 1):
            shell("samtools merge --threads {threads} {params} {output} {input}")
        else:
            shell("ln -s {input} {output}")

rule index_merged_bam:
    input:
        rules.bam_merge.output
    output:
        protected("{assayID}/{runID}/{outdir}/{reference_version}/{application}/{duplicates}/{sample}.bam.bai")
    shell:
        "samtools index {input} {output}"
