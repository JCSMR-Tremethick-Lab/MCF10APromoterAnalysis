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
        trimmed_read1 = "{processed_dir}/{trim_data}/{unit}.end1.fastq.gz",
        trimmed_read2 = "{processed_dir}/{trim_data}/{unit}.end2.fastq.gz",
        report_html = "{processed_dir}/{trim_data}/{unit}_report.html",
        report_json = "{processed_dir}/{trim_data}/{unit}_report.json"
    shell:
        """
            fastp -i {read1} -I {read2} -o {trimmed_read1} -O {trimmed_read2} --html {report_html} --json {report_json} --thread {threads}
        """
