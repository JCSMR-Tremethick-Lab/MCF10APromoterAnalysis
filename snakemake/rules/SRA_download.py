__author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2016-06-08"

# vim: syntax=python tabstop=4 expandtab
# coding: utf-8


"""
Rules for downloading and processing SRA data to FASTQ

For use, include in your workflow.
"""

from snakemake.exceptions import MissingInputException
import snakemake.utils

SRA_IDS, = glob_wildcards("SRA/requested/{id}")

def make_sra_urls(wildcards):
    url = []
    url_prefix = "ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra"
    i = wildcards['sra_ids']
    url.append(url_prefix + "/" + i[0:3] + "/" + i[0:6] + "/" + i + "/" + i + ".sra")
    return(url)

rule download_sra:
    params:
        url = make_sra_urls
    output:
        "SRA/sra/{sra_ids}.sra"
    shell:
        """
            wget {params.url} {output}
        """

rule fastq_dump:
    input:
        rules.download_sra.output
    output:
        "SRA/fastq/{sra_ids}_1.fastq.gz",
        "SRA/fastq/{sra_ids}_2.fastq.gz"
    shell:
        """
            fastq-dump --split-files --gzip {input}
        """

rule download:
    input:
        expand("SRA/fastq/{sra_ids}_1.fastq.gz", sra_ids = SRA_IDS),
        expand("SRA/fastq/{sra_ids}_2.fastq.gz", sra_ids = SRA_IDS)
