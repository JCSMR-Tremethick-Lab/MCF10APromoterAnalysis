__author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2018-06-15"

# vim: syntax=python tabstop=4 expandtab
# coding: utf-8


"""
Rules for running MACS2 for ChIP peak on shor fragments BED files.
This workflow does not include Input as control but rather uses a low stringency
setting in order to detect as many sites of potential short fragment enrichment
prior to performing TFBS motif analysis.
Stand alone workflow.
"""

rule macs2_callpeak:
    params:
        macs2_dir = "/home/sebastian/miniconda3/envs/py27/bin",
        name = lambda wildcards: wildcards.smallFragments
    input:
        chip = "/home/sebastian/Data/Collaborations/FSU/PromoterSeqCap/SmallFragments/{smallFragments}.bed"
    output:
        "/home/sebastian/Data/Collaborations/FSU/PromoterSeqCap/SmallFragments/macs2PeakCalling/{smallFragments}"
    shell:
        """
            {params.macs2_dir}/macs2 callpeak -f BED \
                                              -g 40999507\
                                              -t {input.chip}\
                                              -n {params.name}\
                                              --nomodel\
                                              --extsize 125\
                                              --outdir {output}\
                                              --call-summits\
                                              -p 0.1\
                                              --bdg\
                                              --trackline
        """

rule all:
    input:
        expand("/home/sebastian/Data/Collaborations/FSU/PromoterSeqCap/SmallFragments/macs2PeakCalling/{smallFragments}"
                smallFragments = ["TOTALcombined_A_H2AZ_000-125.bed",
        "TOTALcombined_A_Inp_000-125.bed",
        "TOTALcombined_A_TGFb_H2AZ_000-125.bed",
        "TOTALcombined_A_TGFb_Inp_000-125.bed",
        "TOTALcombined_CA1a_H2AZ_000-125.bed",
        "TOTALcombined_CA1a_Inp_000-125.bed",
        "TOTALcombined_shH2AZ_Inp_000-125.bed"])
