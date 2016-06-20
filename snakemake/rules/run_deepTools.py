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
        expand("{processed_dir}/{genome_version}/deepTools/bamCoverage/{sample}.bw",
               processed_dir = config["processed_dir"],
               genome_version = "hg38",
               sample = config["sample"]),
        expand("{processed_dir}/{genome_version}/deepTools/computeMatrix_referencePoint/{region}.{sample}.matrix.gz",
               processed_dir = config["processed_dir"],
               genome_version = "hg38",
               sample = config["sample"],
               region = ["EMT_markers",
                         "EMT_markers_down",
                         "EMT_markers_up",
                         "seqCapTargets_hg38",
                         "seqCap_Targets_EMT_markers",
                         "seqCap_Targets_EMT_markers_up",
                         "seqCap_Targets_EMT_markers_down"]),
        expand("{processed_dir}/{genome_version}/deepTools/computeMatrix_referencePoint/{region}.{sample}.matrix",
               processed_dir = config["processed_dir"],
               genome_version = "hg38",
               sample = config["sample"],
               region = ["EMT_markers",
                         "EMT_markers_down",
                         "EMT_markers_up",
                         "seqCapTargets_hg38",
                         "seqCap_Targets_EMT_markers",
                         "seqCap_Targets_EMT_markers_up",
                         "seqCap_Targets_EMT_markers_down"]),
        expand("{processed_dir}/{genome_version}/deepTools/plotProfile/profile.{region}.{sample}.pdf",
                processed_dir = config["processed_dir"],
                genome_version = "hg38",
                sample = config["sample"],
                region = ["EMT_markers",
                          "EMT_markers_down",
                          "EMT_markers_up",
                          "seqCapTargets_hg38",
                          "seqCap_Targets_EMT_markers",
                          "seqCap_Targets_EMT_markers_up",
                          "seqCap_Targets_EMT_markers_down"]),
        expand("{processed_dir}/{genome_version}/deepTools/plotProfile/profile.{region}.{sample}.data",
               processed_dir = config["processed_dir"],
               genome_version = "hg38",
               sample = config["sample"],
               region = ["EMT_markers",
                         "EMT_markers_down",
                         "EMT_markers_up",
                         "seqCapTargets_hg38",
                         "seqCap_Targets_EMT_markers",
                         "seqCap_Targets_EMT_markers_up",
                         "seqCap_Targets_EMT_markers_down"]),
        expand("{processed_dir}/{genome_version}/deepTools/plotProfile/profile.{region}.{sample}.bed",
               processed_dir = config["processed_dir"],
               genome_version = "hg38",
               sample = config["sample"],
               region = ["EMT_markers",
                         "EMT_markers_down",
                         "EMT_markers_up",
                         "seqCapTargets_hg38",
                         "seqCap_Targets_EMT_markers",
                         "seqCap_Targets_EMT_markers_up",
                         "seqCap_Targets_EMT_markers_down"])

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
    version:
        0.2
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
                                           --numberOfProcessors 8 \
                                           --normalizeUsingRPKM \
                                           --smoothLength 30 \
                                           --centerReads \
                                           --skipNonCoveredRegions
        """
rule computeMatrix_referencePoint:
    version:
        0.3
    params:
        deepTools_dir = config["deepTools_dir"]
    input:
        file = "{processed_dir}/{genome_version}/deepTools/bamCoverage/{sample}.bw",
        region = "{processed_dir}/{genome_version}/deepTools/regionFiles/{region}.bed"
    output:
        matrix_gz = "{processed_dir}/{genome_version}/deepTools/computeMatrix_referencePoint/{region}.{sample}.matrix.gz",
        matrix = "{processed_dir}/{genome_version}/deepTools/computeMatrix_referencePoint/{region}.{sample}.matrix"
    shell:
        """
        {params.deepTools_dir}/computeMatrix reference-point \
                                             --referencePoint TSS \
                                             --regionsFileName {input.region} \
                                             --scoreFileName {input.file} \
                                             --upstream 300 \
                                             --downstream 300 \
                                             --missingDataAsZero \
                                             --skipZeros \
                                             --outFileName {output.matrix_gz} \
                                             --outFileNameMatrix {output.matrix}
        """

rule plotProfile:
    version:
        0.1
    params:
        deepTools_dir = config["deepTools_dir"]
    input:
        "{processed_dir}/{genome_version}/deepTools/computeMatrix_referencePoint/{region}.{sample}.matrix.gz"
    output:
        figure = "{processed_dir}/{genome_version}/deepTools/plotProfile/profile.{region}.{sample}.pdf",
        data = "{processed_dir}/{genome_version}/deepTools/plotProfile/profile.{region}.{sample}.data",
        regions = "{processed_dir}/{genome_version}/deepTools/plotProfile/profile.{region}.{sample}.bed"
    shell:
        """
            {params.deepTools_dir}/plotProfile --matrixFile {input} \
                                               --outFileName {output.figure} \
                                               --outFileNameData {output.data} \
                                               --outFileSortedRegions {output.regions} \
                                               --plotType se
        """


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
