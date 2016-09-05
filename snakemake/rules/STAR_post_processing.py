__author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2016-09-05"

# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

from snakemake.exceptions import MissingInputException
import os

"""
Rules for running deepTools analysis on RNA-Seq data
For usage, include this in your workflow.
"""

rule all:
    input:
        expand("./{assayID}/{runID}/{processed_dir}/{reference_version}/deepTools/results.npz",
               assayID = "RNA-Seq",
               runID = "NB501086_0067_RDomaschenz_JCSMR_RNASeq",
               processed_dir = config["processed_dir"],
               reference_version = config["references"]["hg19"]["version"]),
        expand("./{assayID}/{runID}/{processed_dir}/{reference_version}/deepTools/multiBamSummary/raw_counts.txt",
               assayID = "RNA-Seq",
               runID = "NB501086_0067_RDomaschenz_JCSMR_RNASeq",
               processed_dir = config["processed_dir"],
               reference_version = config["references"]["hg19"]["version"]),
        expand("./{assayID}/{runID}/{processed_dir}/{reference_version}/deepTools/plotCorrelation/heatmap_SpearmanCorr_readCounts.png",
               assayID = "RNA-Seq",
               runID = "NB501086_0067_RDomaschenz_JCSMR_RNASeq",
               processed_dir = config["processed_dir"],
               reference_version = config["references"]["hg19"]["version"]),
        expand("./{assayID}/{runID}/{processed_dir}/{reference_version}/deepTools/plotPCA/PCA_readCounts.png",
               assayID = "RNA-Seq",
               runID = "NB501086_0067_RDomaschenz_JCSMR_RNASeq",
               processed_dir = config["processed_dir"],
               reference_version = config["references"]["hg19"]["version"])

rule multiBamSummary:
    version:
        0.3
    params:
        deepTools_dir = config["deepTools_dir"],
        binSize = 1000
    input:
        expand("./{assayID}/{runID}/{processed_dir}/{reference_version}/STAR/full/{unit}.aligned.bam",
               assayID = "RNA-Seq",
               runID = "NB501086_0067_RDomaschenz_JCSMR_RNASeq",
               units = config["samples"]["RNA-Seq"],
               processed_dir = config["processed_dir"],
               reference_version = config["hg19"]["version"])
    output:
        npz = "./{assayID}/{runID}/{processed_dir}/{reference_version}/deepTools/results.npz",
        raw = "./{assayID}/{runID}/{processed_dir}/{reference_version}/deepTools/multiBamSummary/raw_counts.txt"
    shell:
        """
        {params.deepTools_dir}/multiBamSummary BED-file --BED seqCapTargets_hg38.bed \
                                                        --bamfiles {input} \
                                                        --numberOfProcessors 8 \
                                                        --centerReads \
                                                        --binSize {params.binSize} \
                                                        --outFileName {output.npz} \
                                                        --outRawCounts {output.raw}
        """

rule plotCorrelation_heatmap:
    version:
        0.3
    params:
        deepTools_dir = config["deepTools_dir"]
    input:
        "./{assayID}/{runID}/{processed_dir}/{reference_version}/deepTools/results.npz"
    output:
        "./{assayID}/{runID}/{processed_dir}/{reference_version}/deepTools/plotCorrelation/heatmap_SpearmanCorr_readCounts.png",
        "./{assayID}/{runID}/{processed_dir}/{reference_version}/deepTools/plotCorrelation/heatmap_SpearmanCorr_readCounts.tab"
    shell:
        """
        {params.deepTools_dir}/plotCorrelation -in {input} \
                                               --corMethod spearman \
                                               --skipZeros \
                                               --plotTitle "Spearman Correlation of Read Counts per genomic bin" \
                                               --whatToPlot heatmap \
                                               --colorMap RdYlBu \
                                               --plotNumbers \
                                               -o {output[0]} \
                                               --outFileCorMatrix {output[1]}
        """

rule plotPCA:
    version:
        0.3
    params:
        deepTools_dir = config["deepTools_dir"]
    input:
        "./{assayID}/{runID}/{processed_dir}/{reference_version}/deepTools/results.npz"
    output:
        "./{assayID}/{runID}/{processed_dir}/{reference_version}/deepTools/plotPCA/PCA_readCounts.png"
    shell:
        """
        {params.deepTools_dir}/plotPCA -in {input} \
                                       -o {output} \
                                       -T "PCA of read counts per genomic bin" \
        """
