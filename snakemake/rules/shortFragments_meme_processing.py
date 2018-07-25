__author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2018-06-15"

# vim: syntax=python tabstop=4 expandtab
# coding: utf-8


"""
Rules for running MEME
"""
import os
home = os.environ['HOME']

rule get_summit_sequences:
    version:
        "1.0"
    params:
        summitsSeqWidth = 500, # recommended by MEME
        peaksMinPileUp = 4,
        peaksMinQval = 2,
        BSgenome = "BSgenome.Hsapiens.UCSC.hg19"
    input:
        summits = "/home/sebastian/Data/Collaborations/FSU/PromoterSeqCap/SmallFragments/macs2PeakCalling/{smallFragments}/{smallFragments}_summits.bed",
        peaks = "/home/sebastian/Data/Collaborations/FSU/PromoterSeqCap/SmallFragments/macs2PeakCalling/{smallFragments}/{smallFragments}_peaks.xls"
    output:
        summitsSeqFile = "/home/sebastian/Data/Collaborations/FSU/PromoterSeqCap/SmallFragments/summitSequences/{smallFragments}_summits.fasta"
    script:
        home + "/Development/JCSMR-Tremethick-Lab/Breast/snakemake/scripts/prepare_summit_sequences.R"


rule run_meme:
    version:
        "1.0"
    params:
        meme_bin = home + "/meme/bin/meme",
        minw = 8,
        maxw = 15,
        nmotifs = 1000,
        evt = 0.05, # e-value threshold
        mod = "zoops"
    threads:
        10
    input:
        summitsSeqFile = home + "/Data/Collaborations/FSU/PromoterSeqCap/SmallFragments/summitSequences/{smallFragments}_summits.fasta",
        promotersHMM = home + "/Data/Collaborations/FSU/PromoterSeqCap/SmallFragments/promoterSequences.hmm"
    output:
        meme_out = home + "/Data/Collaborations/FSU/PromoterSeqCap/SmallFragments/meme/{memeObjectiveFunction}/{smallFragments}",
        memeOutput = home + "/Data/Collaborations/FSU/PromoterSeqCap/SmallFragments/meme/{memeObjectiveFunction}/{smallFragments}/meme.html"
    shell:
        """
            {params.meme_bin} -oc {output.meme_out}\
                              -dna \
                              -bfile {input.promotersHMM} \
                              -minw {params.minw} \
                              -maxw {params.maxw} \
                              -nmotifs {params.nmotifs} \
                              -evt {params.evt} \
                              -p {threads} \
                              -objfun {wildcards.memeObjectiveFunction} \
                              -mod {params.mod} \
                              {input.summitsSeqFile}
        """

rule run_tomtom:
    version:
        "1.0"
    params:
        tomtom = home + "/meme/bin/tomtom"
    threads:
        1
    input:
        motifDB = home + "/Data/References/MEME/motif_databases/HUMAN/HOCOMOCOv11_full_HUMAN_mono_meme_format.meme",
        memeOutput = home + "/Data/Collaborations/FSU/PromoterSeqCap/SmallFragments/meme/{memeObjectiveFunction}/{smallFragments}/meme.html"
    output:
        tomtom_out = home + "/Data/Collaborations/FSU/PromoterSeqCap/SmallFragments/tomtom/{memeObjectiveFunction}/{smallFragments}"
    shell:
        """
            {params.tomtom} -oc {output.tomtom_out} {input.memeOutput} {input.motifDB}
        """

rule motif_summary:
    version:
        "1.0"
    params:
    threads:
        1
    input:
        home + "/Data/Collaborations/FSU/PromoterSeqCap/SmallFragments/meme/{memeObjectiveFunction}/{smallFragments}"
    output:
        home + "/Data/Collaborations/FSU/PromoterSeqCap/SmallFragments/meme/{memeObjectiveFunction}/{smallFragments}/meme.summary.txt"
    shell:
        """
            grep "MOTIF" {input}/meme.txt > {output}
        """



rule all:
    input:
        expand("/home/sebastian/Data/Collaborations/FSU/PromoterSeqCap/SmallFragments/tomtom/{memeObjectiveFunction}/{smallFragments}/meme.summary.txt",
                memeObjectiveFunction = ["cd", "ce"],
                smallFragments = ["TOTALcombined_A_H2AZ_000-125",
                                  "TOTALcombined_A_Inp_000-125",
                                  "TOTALcombined_A_TGFb_H2AZ_000-125",
                                  "TOTALcombined_A_TGFb_Inp_000-125",
                                  "TOTALcombined_CA1a_H2AZ_000-125",
                                  "TOTALcombined_CA1a_Inp_000-125",
                                  "TOTALcombined_shH2AZ_Inp_000-125"])
