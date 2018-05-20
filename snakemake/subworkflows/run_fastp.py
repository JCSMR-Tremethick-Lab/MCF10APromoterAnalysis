__author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2018-05-18"

# vim: syntax=python tabstop=4 expandtab
# coding: utf-8


"""
Rules for trimming reads with fastq
(https://github.com/OpenGene/fastp)

For usage, include this in your workflow.
"""

rule run_fastp:
    version:
        0.1
    threads:
        4
    input:
        read1 = "fastq/{unit}.end1.fastq.gz",
        read2 = "fastq/{unit}.end2.fastq.gz"
    output:
        trimmed_read1 = "fastq/{unit}.end1.trimmed.fastq.gz",
        trimmed_read2 = "fastq/{unit}.end2.trimmed.fastq.gz",
        report_html = "fastq/{unit}_report.html",
        report_json = "fastq/{unit}_report.json"
    shell:
        "fastp -i {input.read1} -I {input.read2} -o {output.trimmed_read1} -O {output.trimmed_read2} --html {output.report_html} --json {output.report_json} --thread {threads}"

rule bowtie2_pe:
    version:
        "0.2"
    params:
        max_in = 250,
        bt2_index = home + config["references"][REF_GENOME]["bowtie2"][REF_VERSION]
    threads:
        8
    input:
        trimmed_read1 = "fastq/{unit}.end1.trimmed.fastq.gz",
        trimmed_read2 = "fastq/{unit}.end2.trimmed.fastq.gz"
    output:
        "{outdir}/{reference_version}/bowtie2/{unit}.bam"
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
