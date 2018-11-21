__author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2018-06-15"

# vim: syntax=python tabstop=4 expandtab
# coding: utf-8
from snakemake.utils import min_version
min_version("5.0")

"""
Rules for running MEME
"""
import os
home = os.environ['HOME']


rule run_meme:
    version:
        "1.0"
    params:
        meme_bin = home + "/meme/bin/meme",
        minw = 8,
        maxw = 15,
        nmotifs = 1000,
        evt = 0.05, # e-value threshold
        mod = "zoops",
        nbrief = 2000
    threads:
        10
    input:
        summitsSeqFile = home + "/Data/Collaborations/FSU/PromoterSeqCap/SmallFragments/summitSequences/{smallFragments}_summits.fasta",
        promotersHMM = home + "/Data/Collaborations/FSU/PromoterSeqCap/SmallFragments/promoterSequences.hmm"
    output:
        directory(home + "/Data/Collaborations/FSU/PromoterSeqCap/SmallFragments/meme/{memeObjectiveFunction}/{smallFragments}")
#        memeOutput = home + "/Data/Collaborations/FSU/PromoterSeqCap/SmallFragments/meme/{memeObjectiveFunction}/{smallFragments}/meme.html"
    shell:
        """
            {params.meme_bin} -oc {output[0]}\
                              -dna \
                              -bfile {input.promotersHMM} \
                              -minw {params.minw} \
                              -maxw {params.maxw} \
                              -nmotifs {params.nmotifs} \
                              -evt {params.evt} \
                              -p {threads} \
                              -objfun {wildcards.memeObjectiveFunction} \
                              -mod {params.mod} \
                              -brief {params.nbrief} \
                              {input.summitsSeqFile}
        """


rule run_tomtom:
    version:
        "1.0"
    params:
        tomtom_bin = home + "/meme/bin/tomtom"
    threads:
        1
    input:
        motifDB = home + "/Data/References/MEME/motif_databases/HUMAN/HOCOMOCOv11_full_HUMAN_mono_meme_format.meme",
        memeOutput = home + "/Data/Collaborations/FSU/PromoterSeqCap/SmallFragments/meme/{memeObjectiveFunction}/{smallFragments}/meme.html"
    output:
        directory(home + "/Data/Collaborations/FSU/PromoterSeqCap/SmallFragments/tomtom/{memeObjectiveFunction}/{smallFragments}")
    shell:
        """
            {params.tomtom_bin} -oc {output[0]} {input.memeOutput} {motifDB}
        """


rule motif_summary:
    version:
        "1.0"
    threads:
        1
    input:
        home + "/Data/Collaborations/FSU/PromoterSeqCap/SmallFragments/meme/{memeObjectiveFunction}/{smallFragments}"
    output:
        home + "/Data/Collaborations/FSU/PromoterSeqCap/SmallFragments/meme_summary/{memeObjectiveFunction}/{smallFragments}/summary.txt"
    shell:
        """
            grep "MOTIF" {input}/meme.txt > {output}
        """


rule run_fimo: # to get sequences in which motifs can be found
    version:
        "1.0"
    params:
        bfile = home + "/Data/Collaborations/FSU/PromoterSeqCap/SmallFragments/promoterSequences.hmm"
    threads:
        1
    input:
        meme = home + "/Data/Collaborations/FSU/PromoterSeqCap/SmallFragments/meme/{memeObjectiveFunction}/{smallFragments}/meme.html",
        fasta = home + "/Data/Collaborations/FSU/PromoterSeqCap/SmallFragments/summitSequences/{smallFragments}_summits.fasta"
    output:
        directory(home + "/Data/Collaborations/FSU/PromoterSeqCap/SmallFragments/fimo/{memeObjectiveFunction}/{smallFragments}")
    shell:
        """
             ~/meme/bin/fimo -bfile {params.bfile} \
                             -oc {output} \
                             {input.meme} \
                             {input.fasta}
        """


rule all:
    input:
        expand("/home/sebastian/Data/Collaborations/FSU/PromoterSeqCap/SmallFragments/fimo/{memeObjectiveFunction}/{smallFragments}",
                memeObjectiveFunction = ["cd", "ce"],
                smallFragments = ["TOTALcombined_A_H2AZ_000-125",
                                  "TOTALcombined_A_Inp_000-125",
                                  "TOTALcombined_A_TGFb_H2AZ_000-125",
                                  "TOTALcombined_A_TGFb_Inp_000-125",
                                  "TOTALcombined_CA1a_H2AZ_000-125",
                                  "TOTALcombined_CA1a_Inp_000-125",
                                  "TOTALcombined_shH2AZ_Inp_000-125"])