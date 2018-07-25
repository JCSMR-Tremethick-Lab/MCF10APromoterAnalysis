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
import os
home = os.environ['HOME']

rule macs2_callpeak:
    params:
        macs2_dir = "/home/sebastian/miniconda3/envs/py27/bin",
        name = lambda wildcards: wildcards.smallFragments
    input:
        chip = "/home/sebastian/Data/Collaborations/FSU/PromoterSeqCap/SmallFragments/{smallFragments}.bed"
    output:
        "/home/sebastian/Data/Collaborations/FSU/PromoterSeqCap/SmallFragments/macs2PeakCalling/{smallFragments}",
        "/home/sebastian/Data/Collaborations/FSU/PromoterSeqCap/SmallFragments/macs2PeakCalling/{smallFragments}_summits.bed",
        "/home/sebastian/Data/Collaborations/FSU/PromoterSeqCap/SmallFragments/macs2PeakCalling/{smallFragments}_peaks.xls"
    shell:
        """
            {params.macs2_dir}/macs2 callpeak -f BED \
                                              -g 40999507\
                                              -t {input.chip}\
                                              -n {params.name}\
                                              --nomodel\
                                              --extsize 125\
                                              --outdir {output[0]}\
                                              --call-summits\
                                              -p 0.1\
                                              --bdg\
                                              --trackline
        """


rule get_summit_sequences:
    version:
        "1.0"
    params:
        summitsSeqWidth = "500", # recommended by MEME
        peaksMinPileUp = "4",
        peaksMinQval = "2",
        BSgenome = "BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19"
    input:
        summits = "/home/sebastian/Data/Collaborations/FSU/PromoterSeqCap/SmallFragments/macs2PeakCalling/{smallFragments}_summits.bed",
        peaks = "/home/sebastian/Data/Collaborations/FSU/PromoterSeqCap/SmallFragments/macs2PeakCalling/{smallFragments}_peaks.xls"
    output:
        summitsSeqFile = "/home/sebastian/Data/Collaborations/FSU/PromoterSeqCap/SmallFragments/macs2PeakCalling/{smallFragments}_summits.fasta"
    script:
        home + "/Development/JCSMR-Tremethick-Lab/Breast/snakemake/scripts/prepare_summit_sequences.R"

rule sort_bedGraph:
    input:
        macs2output = "/home/sebastian/Data/Collaborations/FSU/PromoterSeqCap/SmallFragments/macs2PeakCalling/{smallFragments}",
        bdg = "/home/sebastian/Data/Collaborations/FSU/PromoterSeqCap/SmallFragments/macs2PeakCalling/{smallFragments}/{smallFragments}_treat_pileup.bdg"
    output:
        temp("/home/sebastian/Data/Collaborations/FSU/PromoterSeqCap/SmallFragments/macs2PeakCalling/{smallFragments}/{smallFragments}_treat_pileup.temp")
    shell:
        """
            grep -v "chrM" {input.bdg} | grep -v "track" | sort -k1,1 -k2,2n - > {output}
        """


rule bdg_to_bigWig:
    params:
        macs2_dir = "/home/sebastian/miniconda3/envs/py27/bin",
        name = lambda wildcards: wildcards.smallFragments,
        chromSizes = "~/Data/References/Genomes/Homo_sapiens/GRCh37_hg19_UCSC/chromSizes.txt"
    input:
        macs2output = "/home/sebastian/Data/Collaborations/FSU/PromoterSeqCap/SmallFragments/macs2PeakCalling/{smallFragments}",
        bdg = "/home/sebastian/Data/Collaborations/FSU/PromoterSeqCap/SmallFragments/macs2PeakCalling/{smallFragments}/{smallFragments}_treat_pileup.temp"
    output:
        bw = "/home/sebastian/Data/Collaborations/FSU/PromoterSeqCap/SmallFragments/macs2PeakCalling/{smallFragments}/{smallFragments}_treat_pileup.bw"
    shell:
        """
            bedGraphToBigWig {input.bdg} {params.chromSizes} {output.bw}
        """

rule all:
    input:
        expand("/home/sebastian/Data/Collaborations/FSU/PromoterSeqCap/SmallFragments/macs2PeakCalling/{smallFragments}/{smallFragments}{suffix}",
                smallFragments = ["TOTALcombined_A_H2AZ_000-125",
                                  "TOTALcombined_A_Inp_000-125",
                                  "TOTALcombined_A_TGFb_H2AZ_000-125",
                                  "TOTALcombined_A_TGFb_Inp_000-125",
                                  "TOTALcombined_CA1a_H2AZ_000-125",
                                  "TOTALcombined_CA1a_Inp_000-125",
                                  "TOTALcombined_shH2AZ_Inp_000-125"],
                suffix = ["_treat_pileup.bw", "_summits.fasta"])
