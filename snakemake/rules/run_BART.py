__author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2018-06-26"

# vim: syntax=python tabstop=4 expandtab
# coding: utf-8


"""
Rule for processing of small fragment BED files.

to be used stand alone.
"""

import os
import fnmatch
from snakemake.exceptions import MissingInputException

# set local variables
home = os.environ['HOME']

rule bart_profile:
    threads:
        8
    params:
    input:
        BED="BEDs/{sample}.bed"
    output:
        dir="{sample}"
    run:
        "bart profile -i {input.BED} -f bed -s hg38 -p {threads} --outdir {output}"

rule all:
    input:
        ["TOTALcombined_A_H2AZ_000-100",
        "TOTALcombined_A_H2AZ_000-125",
        "TOTALcombined_A_Inp_000-100",
        "TOTALcombined_A_Inp_000-125",
        "TOTALcombined_A_TGFb_H2AZ_000-100",
        "TOTALcombined_A_TGFb_H2AZ_000-125",
        "TOTALcombined_A_TGFb_Inp_000-100",
        "TOTALcombined_A_TGFb_Inp_000-125",
        "TOTALcombined_CA1a_H2AZ_000-100",
        "TOTALcombined_CA1a_H2AZ_000-125",
        "TOTALcombined_CA1a_Inp_000-100",
        "TOTALcombined_CA1a_Inp_000-125",
        "TOTALcombined_shH2AZ_Inp_000-100",
        "TOTALcombined_shH2AZ_Inp_000-125"]
