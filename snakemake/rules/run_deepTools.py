__author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2016-03-01"

# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

from snakemake.exceptions import MissingInputException
import os

"""
Rules for running deepTools analysis on ChIP-Seq data
For usage, include this in your workflow.
"""

rule all:
    input:
        "processed_data/hg38/deepTools/results.npz",
        # expand("deepTools/bamPEFragmentSize/{samples}_histogram.png", samples = config["units"]),
        "processed_data/hg38/deepTools/plotCorrelation/heatmap_SpearmanCorr_readCounts.png",
        "processed_data/hg38/deepTools/plotPCA/PCA_readCounts.png",
        expand("{processed_dir}/{genome_version}/deepTools/bamCoverage/{sample}.bw", processed_dir = config["processed_dir"], genome_version = "hg38", sample = config["sample"])
        # expand("deepTools/plotFingerprint/{dup}_fingerprints.{dup_suff}.png", dup = "duplicates_removed", dup_suff = "DeDup"),
        # expand("./deepTools/bamCompare/{chip}_vs_Input.{norm}.bw", chip = ("H2AZ", "H2ABbd"), norm = ("SES", "readCount")),
        # expand("deepTools/computeMatrix_referencePoint/{region}.{sample}.{norm}.{type}", region = ("ctaGenes", "allGenes", "ctaGenesExpressed"), sample = ("H2ABbd_vs_Input", "H2AZ_vs_Input"), norm = ("SES", "readCount"), type = ("matrix.gz", "matrix")),
        # expand("deepTools/computeMatrix_scaleRegions/{region}.{sample}.{norm}.{type}", region = ("ctaGenes", "allGenes", "ctaGenesExpressed"), sample = ("H2ABbd_vs_Input", "H2AZ_vs_Input"), norm = ("SES", "readCount"), type = ("matrix.gz", "matrix")),
        # expand("deepTools/plotHeatmap/{matrix_dir}/heatmap.{region}.{sample}.{norm}.{type}", region = ("ctaGenes", "allGenes", "ctaGenesExpressed"), matrix_dir = ("computeMatrix_referencePoint"), sample = ("H2ABbd_vs_Input", "H2AZ_vs_Input"), norm = ("SES", "readCount"), type = ("pdf", "data", "bed", "matrix"))

rule multiBamSummary:
    params:
        deepTools_dir = config["deepTools_dir"]
    input:
        expand("./processed_data/duplicates_removed/{units}.DeDup.sorted.fastq_q20.bam", units = config["units"])
    output:
        "deepTools/results.npz"
    shell:
        """
        {params.deepTools_dir}/multiBamSummary BED-file --BED seqCapTargets_hg38.bed \
                                                        --bamfiles {input} \
                                                        --numberOfProcessors max \
                                                        --centerReads \
                                                        --outFileName {output}
        """

rule plotCorrelation_heatmap:
    params:
        deepTools_dir = config["deepTools_dir"]
    input:
        "deepTools/results.npz"
    output:
        "deepTools/plotCorrelation/heatmap_SpearmanCorr_readCounts.png",
        "deepTools/plotCorrelation/heatmap_SpearmanCorr_readCounts.tab"
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
    params:
        deepTools_dir = config["deepTools_dir"]
    input:
        "deepTools/results.npz"
    output:
        "deepTools/plotPCA/PCA_readCounts.png"
    shell:
        """
        {params.deepTools_dir}/plotPCA -in {input} \
                                       -o {output} \
                                       -T "PCA of read counts per genomic bin" \
        """

rule bamPEFragmentSize:
    params:
        deepTools_dir = config["deepTools_dir"]
    input:
        "./processed_data/duplicates_marked/{units}.Q20.sorted.MkDup.bam"
    output:
        "deepTools/bamPEFragmentSize/{units}_histogram.png"
    shell:
        """
        {params.deepTools_dir}/bamPEFragmentSize --histogram {output} {input}
        """

rule plotFingerprint:
    params:
        deepTools_dir = config["deepTools_dir"]
    input:
        expand("./processed_data/duplicates_removed/{units}.DeDup.sorted.fastq_q20.bam", units = config["units"])
    output:
        "deepTools/plotFingerprint/duplicates_removed_fingerprints.png"
    shell:
        """
        {params.deepTools_dir}/plotFingerprint --bamfiles {input} \
                                               --numberOfProcessors max \
                                               --centerReads \
                                               --plotTitle "Library complexity" \
                                               --skipZeros \
                                               --plotFile {output}
        """

rule bamCoverage_MNase:
    params:
        deepTools_dir = config["deepTools_dir"]
    input:
        "{processed_dir}/{genome_version}/duplicates_removed/{sample}.Q20.DeDup.sorted.bam"
    output:
        "{processed_dir}/{genome_version}/deepTools/bamCoverage/{sample}.bw"
    shell:
        """
        {params.deepTools_dir}/bamCoverage --bam {input} \
                                           --outFileName {output} \
                                           --outFileFormat bigwig \
                                           --MNase \
                                           --binSize 10 \
                                           --numberOfProcessors 4 \
                                           --normalizeUsingRPKM \
                                           --smoothLength 20 \
                                           --centerReads \
                                           --skipNonCoveredRegions
        """

# rule bamCompare:
#     params:
#         deepTools_dir = config["deepTools_dir"],
#     input:
#         control = expand("./processed_data/{dup}/{samples}.Q20.sorted.{dup_suff}.bam", samples = "Input", dup = "duplicates_removed", dup_suff = "DeDup"),
#         chip = "./processed_data/duplicates_removed/{chip}.Q20.sorted.DeDup.bam"
#     output:
#         file = "./deepTools/bamCompare/{chip}_vs_Input.{norm}.bw",
#     shell:
#         """
#         {params.deepTools_dir}/bamCompare --bamfile1 {input.chip} \
#                                           --bamfile2 {input.control} \
#                                           --outFileName {output.file} \
#                                           --outFileFormat bigwig \
#                                           --scaleFactorsMethod {wildcards.norm} \
#                                           --ratio log2 \
#                                           --numberOfProcessors max \
#                                           --skipNonCoveredRegions
#         """
#
# rule computeMatrix_scaleRegions:
#     version:
#         0.2
#     params:
#         deepTools_dir = config["deepTools_dir"]
#     input:
#         files = expand("deepTools/{data_dir}/{samples}.bw", samples = config["units"], data_dir = "bamCoverage"),
#         regions = "deepTools/regionFiles/{region}.bed"
#     output:
#         matrix_gz = "deepTools/computeMatrix_scaleRegions/{region}.{samples}.{norm}.matrix.gz",
#         matrix = "deepTools/computeMatrix_scaleRegions/{region}.{samples}.{norm}.matrix"
#     shell:
#         """
#         {params.deepTools_dir}/computeMatrix scale-regions \
#                                              --regionsFileName {input.regions} \
#                                              --scoreFileName {input.files} \
#                                              --regionBodyLength 5000 \
#                                              --upstream 1500 \
#                                              --downstream 1500 \
#                                              --unscaled5prime 200 \
#                                              --unscaled3prime 200 \
#                                              --missingDataAsZero \
#                                              --outFileName {output.matrix_gz} \
#                                              --outFileNameMatrix {output.matrix}
#         """
#
# rule computeMatrix_referencePoint:
#     version:
#         0.2
#     params:
#         deepTools_dir = config["deepTools_dir"]
#     input:
#         files = "deepTools/bamCompare/{samples}.{norm}.bw",
#         regions = "deepTools/regionFiles/{region}.bed"
#     output:
#         matrix_gz = "deepTools/computeMatrix_referencePoint/{region}.{samples}.{norm}.matrix.gz",
#         matrix = "deepTools/computeMatrix_referencePoint/{region}.{samples}.{norm}.matrix"
#     shell:
#         """
#         {params.deepTools_dir}/computeMatrix reference-point \
#                                              --referencePoint TSS \
#                                              --regionsFileName {input.regions} \
#                                              --scoreFileName {input.files} \
#                                              --upstream 1500 \
#                                              --downstream 1500 \
#                                              --missingDataAsZero \
#                                              --skipZeros \
#                                              --outFileName {output.matrix_gz} \
#                                              --outFileNameMatrix {output.matrix}
#         """
#
# rule plotHeatmap:
#     version:
#         0.2
#     params:
#         deepTools_dir = config["deepTools_dir"]
#     input:
#         "deepTools/{matrix_dir}/{region}.{sample}.{norm}.matrix.gz"
#     output:
#         figure = "deepTools/plotHeatmap/{matrix_dir}/heatmap.{region}.{sample}.{norm}.pdf",
#         data = "deepTools/plotHeatmap/{matrix_dir}/heatmap.{region}.{sample}.{norm}.data",
#         regions = "deepTools/plotHeatmap/{matrix_dir}/heatmap.{region}.{sample}.{norm}.bed",
#         matrix = "deepTools/plotHeatmap/{matrix_dir}/heatmap.{region}.{sample}.{norm}.matrix"
#     shell:
#         """
#         {params.deepTools_dir}/plotHeatmap --matrixFile {input} \
#                                            --outFileName {output.figure} \
#                                            --outFileNameData {output.data} \
#                                            --outFileSortedRegions {output.regions} \
#                                            --outFileNameMatrix {output.matrix} \
#                                            --kmeans 4
#         """
