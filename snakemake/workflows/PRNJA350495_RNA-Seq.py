import pandas as pd
import os
from snakemake.exceptions import MissingInputException
from snakemake.utils import validate, min_version
# set minimum snakemake version #
min_version("5.1.2")

__author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2018-05-22"

rule:
    version: 0.1

localrules:
    all

wrapper_dir = os.environ['HOME'] + "/Development/snakemake-wrappers/bio"
include_prefix= os.environ['HOME'] + "/Development/JCSMR-Tremethick-Lab/Breast/snakemake/rules/"
config_dir = os.environ['HOME'] + "/Development/JCSMR-Tremethick-Lab/Breast/snakemake/configs/"
home = os.environ['HOME']

configfile: config_dir + "config_RNA-Seq.json"
# validate(config, schema="schemas/config.schema.yaml")

samples = pd.read_table(config_dir + "PRNJA350495_samples.tsv").set_index("sample", drop=False)
#validate(samples, schema="schemas/samples.schema.yaml")

units = pd.read_table(config_dir + "PRNJA350495_units.tsv", dtype=str).set_index(["sample", "unit"], drop=False)
units.index = units.index.set_levels([i.astype(str) for i in units.index.levels])  # enforce str in index
#validate(units, schema="schemas/units.schema.yaml")

# functions #
def get_fastq(wildcards):
    return units.loc[(wildcards["sample"]), ["fq1", "fq2"]].dropna()

# targets #
rule all:
    input:
        expand("{outdir}/{reference_version}/kallisto/{sample}",
               outdir = config["processed_dir"],
               reference_version = "GRCh37_hg19_UCSC",
               sample = samples["sample"].tolist())

# set up report #
report: "report/workflow.rst"

rule fastp_pe:
    version:
        0.1
    threads:
        4
    input:
        read1 = lambda wildcards: config["raw_dir"] + "/" + "".join(get_fastq(wildcards)["fq1"].tolist()),
        read2 = lambda wildcards: config["raw_dir"] + "/" + "".join(get_fastq(wildcards)["fq2"].tolist())
    output:
        trimmed_read1 = "trimmed/{sample}.end1.fastq.gz",
        trimmed_read2 = "trimmed/{sample}.end2.fastq.gz",
        report_html = "trimmed/{sample}_report.html",
        report_json = "trimmed/{sample}_report.json"
    shell:
        "fastp -i {input.read1} -I {input.read2} -o {output.trimmed_read1} -O {output.trimmed_read2} --html {output.report_html} --json {output.report_json} --thread {threads}"

rule kallisto_quant:
    message:
        "Running kallisto quant..."
    threads:
        4
    params:
        bootstraps = config["kallisto"]["bootstraps"],
        trim_dir = config["trim_dir"]
    input:
        read1 = "trimmed/{sample}.end1.fastq.gz",
        read2 = "trimmed/{sample}.end1.fastq.gz",
        ki = lambda wildcards: config["references"]["hg19"]["kallisto"][wildcards.reference_version]
    output:
        protected("{processed_dir}/{reference_version}/kallisto/{sample}")
    shell:
        """
            kallisto quant --index={input.ki} \
                           --output-dir={output} \
                           --threads={threads} \
                           --bootstrap-samples={params.bootstraps} \
                           {input.read1} {input.read2}
        """
