__author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2016-12-05"

import os
import fnmatch
from snakemake.exceptions import MissingInputException

rule:
    version: 0.2

localrules:
    all

home = os.environ['HOME']

wrapper_dir = home + "/Development/snakemake-wrappers/bio"

include_prefix = home + "/Development/JCSMR-Tremethick-Lab/Breast/snakemake/rules/"
subworkflow_prefix = home + "/Development/JCSMR-Tremethick-Lab/Breast/snakemake/subworkflows/"

#include:
#    include_prefix + "run_fastp.py"

# define global variables such as reference version of genome so that it can be accessed
# throughout the whole worfklow
REF_GENOME = config["references"]["genomes"][1]

# set local variables
home = os.environ['HOME']
REF_GENOME = "hg19"
REF_VERSION = "GRCh37_hg19_UCSC"

subworkflow read_trimming:
    workdir: "."
    snakefile: subworkflow_prefix + "run_fastp.py"

rule bowtie2_pe:
    version:
        "0.2"
    params:
        max_in = 250,
        bt2_index = home + config["references"][REF_GENOME]["bowtie2"][REF_VERSION]
    threads:
        8
    input:
        trimmed_read1 = read_trimming("fastq/{unit}.end1.trimmed.fastq.gz"),
        trimmed_read2 = read_trimming("fastq/{unit}.end2.trimmed.fastq.gz")
    output:
        temp("{outdir}/{reference_version}/bowtie2/{unit}.bam")
    shell:
        """
            bowtie2 \
            -x {params.bt2_index}\
            --no-mixed \
            --no-discordant \
            --maxins {params.max_in} \
            --threads {threads}\
            --rg-id '{wildcards.unit}' \
            --rg 'LB:{wildcards.unit}' \
            --rg 'SM:{wildcards.unit}' \
            --rg 'PL:Illumina' \
            --rg 'PU:NA' \
            -1 {input.trimmed_read1} \
            -2 {input.trimmed_read2} \
            | samtools view -Sb - > {output}
        """

rule bam_quality_filter:
    params:
        qual = config["alignment_quality"]
    input:
        rules.bowtie2_pe.output
    output:
        temp("{outdir}/{reference_version}/bowtie2/{unit}.qfiltered.bam")
    shell:
        "samtools view -b -h -q {params.qual} {input} > {output}"

rule bam_sort:
    params:
        qual = config["alignment_quality"]
    threads:
        4
    input:
        rules.bam_quality_filter.output
    output:
        temp("{outdir}/{reference_version}/bowtie2/{unit}.qfiltered.sorted.bam")
    shell:
        "samtools sort -@ {threads} {input} -T {wildcards.unit}.Q{params.qual}.sorted -o {output}"

rule bam_mark_duplicates:
    params:
        qual = config["alignment_quality"],
        picard = home + config["picard"],
        temp = home + config["temp_dir"]
    input:
        rules.bam_sort.output
    output:
        temp("{outdir}/{reference_version}/bowtie2/{unit}.qfiltered.sorted.duplicates_marked.bam")
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

rule bam_rmdup:
    input:
        rules.bam_mark_duplicates.output
    output:
        protected("{outdir}/{reference_version}/bowtie2/{unit}.final.bam")
    shell:
        "samtools rmdup {input} {output}"

rule bam_rmdup_index:
    params:
        qual = config["alignment_quality"]
    input:
        rules.bam_rmdup.output
    output:
        protected("{outdir}/{reference_version}/bowtie2/{unit}.final.bam.bai")
    shell:
        "samtools index {input} {output}"

rule run_bt2:
    input:
        expand("{outdir}/{reference_version}/bowtie2/{unit}.bam",
                outdir = config["processed_dir"],
                reference_version = config["references"][REF_GENOME]["version"],
                unit = config["units"])

rule all:
    input:
        expand("{outdir}/{reference_version}/bowtie2/{unit}.final.{suffix}",
               outdir = config["processed_dir"],
               reference_version = config["references"][REF_GENOME]["version"],
               duplicates = "duplicates_removed",
               unit = config["units"],
               qual = config["alignment_quality"],
               suffix = ["bam", "bam.bai"])
