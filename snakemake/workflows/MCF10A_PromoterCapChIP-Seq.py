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
home = os.environ['HOME']
REF_GENOME = "hg19"
REF_VERSION = "GRCh37_hg19_UCSC"

subworkflow read_trimming:
    workdir: "."
    snakefile: subworkflow_prefix + "run_fastp.py"

rule bam_quality_filter:
    params:
        qual = config["alignment_quality"]
    input:
        read_trimming("{outdir}/{reference_version}/bowtie2/{unit}.bam")
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

rule bam_merge:
    version:
        0.4
    params:
        cwd = os.getcwd()
    threads:
        8
    input:
        bam = lambda wildcards: expand("{outdir}/{reference_version}/bowtie2/{unit}.final.{suffix}",
                                 outdir = wildcards["outdir"],
                                 reference_version = wildcards["reference_version"],
                                 unit = config["samples"][wildcards["condition"]][wildcards["type"]][wildcards["sample"] + "_" + wildcards["type"]],
                                 suffix = ["bam"]),
        index = lambda wildcards: expand("{outdir}/{reference_version}/bowtie2/{unit}.final.{suffix}",
                                 outdir = wildcards["outdir"],
                                 reference_version = wildcards["reference_version"],
                                 unit = config["samples"][wildcards["condition"]][wildcards["type"]][wildcards["sample"] + "_" + wildcards["type"]],
                                 suffix = ["bam.bai"])
    output:
        protected("{outdir}/{reference_version}/bowtie2/merged/{sample}_{type}.{condition}.bam")
    run:
        if (len(input.bam) > 1):
            shell("samtools merge --threads {threads} {output} {input.bam}")
        else:
            shell("ln -s {params.cwd}/{input.bam} {output}")

rule index_merged:
    version:
        0.1
    input:
        rules.bam_merge.output
    output:
        protected("{outdir}/{reference_version}/bowtie2/merged/{sample}_{type}.{condition}.bam.bai")
    shell:
        "samtools index {input} {output}"

rule bamCoverage:
    version:
        0.2
    params:
        deepTools_dir = home + config["deepTools_dir"],
        ignore = config["program_parameters"]["deepTools"]["ignoreForNormalization"],
        program_parameters = cli_parameters_bamCoverage
    threads:
        8
    input:
        bam = "{outdir}/{reference_version}/bowtie2/merged/{sample}_{type}.{condition}.bam",
        index = "{outdir}/{reference_version}/bowtie2/merged/{sample}_{type}.{condition}.bam.bai"
    output:
        protected("{outdir}/{reference_version}/{application}/{tool}/{mode}/{normalization}/{sample}_{type}.{condition}.bw")
    shell:
        """
        {params.deepTools_dir}/bamCoverage --bam {input.bam} \
                                           --outFileName {output} \
                                           --outFileFormat bigwig \
                                           {params.program_parameters} \
                                           --numberOfProcessors {threads} \
                                           --normalizeUsingRPKM \
                                           --ignoreForNormalization {params.ignore}
        """

rule all:
    input:
        expand("{outdir}/{reference_version}/{application}/{tool}/{mode}/{normalization}/{sample}_{type}.{condition}.bw",
               outdir = config["processed_dir"],
               reference_version = "GRCh37_hg19_UCSC",
               application = "deepTools",
               tool = "bamCoverage",
               mode = "MNase",
               normalization = "RPKM",
               sample = config["samples"]["sample"],
               type = "Input",
               condition = "Total",
               suffix = ["bam", "bam.bai"]),
        expand("{outdir}/{reference_version}/{application}/{tool}/{mode}/{normalization}/{sample}_{type}.{condition}.bw",
               outdir = config["processed_dir"],
               reference_version = "GRCh37_hg19_UCSC",
               application = "deepTools",
               tool = "bamCoverage",
               mode = "MNase",
               normalization = "RPKM",
               sample = ["MCF10A_WT", "MCF10A_TGFb", "MCF10CA1a_WT"],
               type = "H2AZ",
               condition = "Total",
               suffix = ["bam", "bam.bai"]),
