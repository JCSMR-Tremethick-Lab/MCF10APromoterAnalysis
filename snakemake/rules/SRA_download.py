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
import wget


def make_sra_urls(wildcards):
    url_prefix = "ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra"
    url = []
    i = wildcards["sra_ids"]:
    print(str(i))
    url.append(url_prefix + "/" + i[0:3] + "/" + i[0:6] + "/" + i + "/" + i + ".sra")

rule download_sra:
    input:
        make_sra_urls
    output:
        "SRA/sra/{sra_ids}.sra"
    shell:
        """
            wget {input} --output-file {output}
        """

rule download:
    input:
        expand("SRA/sra/{sra_ids}.sra", sra_ids = ["SRR786494", "SRR786495"])
