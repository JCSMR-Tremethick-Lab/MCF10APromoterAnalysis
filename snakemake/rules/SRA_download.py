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
    url_prefix = "ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra"
    url = []
    for i in wildcards["sra_ids"]:
        url.append(url_prefix + "/" + i[0:3] + "/" + i[0:6] + "/" + i + "/" + i + ".sra")
    return(url)

rule download_sra:
    input:
        "SRA/requested/{sra_ids}"
    output:
        "SRA/sra/{sra_ids}.sra"
    shell:
        """
            prefetch {wildcards.sra_ids} --ascp-options --move-after-transfer=./SRA/sra/
        """

rule download:
    input:
        expand("SRA/sra/{sra_ids}.sra", sra_ids = SRA_IDS)
