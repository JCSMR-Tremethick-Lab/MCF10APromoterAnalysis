__author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2017-01-11"

# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

from snakemake.exceptions import MissingInputException
import os

"""
Rules for running deepTools QC/QC on ChIP-Seq data
For usage, include this in your workflow.
"""

# global functions
def get_sample_labels(wildcards):
    sl = []
    runIDs = config["samples"][wildcards.assayID].keys()
    for i in runIDs:
        for k in config["samples"][wildcards.assayID][i].keys():
            sl.append(k)
    return(sl)


rule multiBamSummary:
    version:
        0.2
    params:
        deepTools_dir = home + config["deepTools_dir"],
        binSize = config["program_parameters"]["deepTools"]["binSize"],
        labels = get_sample_labels
    threads:
        24
    input:
        expand("{assayID}/{runID}/{outdir}/{reference_version}/bowtie2/duplicates_marked/{unit}.Q{qual}.sorted.MkDup.{suffix}",
               assayID = "ChIP-Seq",
               runID = "NB501086_0086_DSTremethick_JCSMR_MCF10A_ChIPseq",
               outdir = config["processed_dir"],
               reference_version = config["references"][REF_GENOME]["version"][0],
               unit = config["samples"]["ChIP-Seq"]["merged"],
               qual = config["alignment_quality"],
               suffix = ["bam"])
    output:
        npz = "{assayID}/{outdir}/{reference_version}/deepTools/multiBamSummary/duplicates_marked/results.npz"
    shell:
        """
            {params.deepTools_dir}/multiBamSummary bins --bamfiles {input} \
                                                        --labels {params.labels} \
                                                        --numberOfProcessors {threads} \
                                                        --centerReads \
                                                        --binSize {params.binSize} \
                                                        --outFileName {output.npz}
        """

rule multiBamSummary_deduplicated:
    version:
        0.2
    params:
        deepTools_dir = home + config["deepTools_dir"],
        binSize = config["program_parameters"]["deepTools"]["binSize"],
        labels = get_sample_labels
    threads:
        24
    input:
        expand("{assayID}/{runID}/{outdir}/{reference_version}/bowtie2/duplicates_removed/{unit}.Q{qual}.sorted.DeDup.{suffix}",
               assayID = "ChIP-Seq",
               runID = "NB501086_0086_DSTremethick_JCSMR_MCF10A_ChIPseq",
               outdir = config["processed_dir"],
               reference_version = config["references"][REF_GENOME]["version"][0],
               unit = config["samples"]["ChIP-Seq"]["merged"],
               qual = config["alignment_quality"],
               suffix = ["bam"])
    output:
        npz = "{assayID}/{outdir}/{reference_version}/deepTools/multiBamSummary/duplicates_removed/results.npz"
    shell:
        """
            {params.deepTools_dir}/multiBamSummary bins --bamfiles {input} \
                                                        --labels {params.labels} \
                                                        --numberOfProcessors {threads} \
                                                        --centerReads \
                                                        --binSize {params.binSize} \
                                                        --outFileName {output.npz}
        """


rule plotCorrelation_heatmap:
    params:
        deepTools_dir = home + config["deepTools_dir"],
        plotTitle = lambda wildcards: "Correlation heatmap - " + wildcards.duplicates
    input:
        npz = "{assayID}/{outdir}/{reference_version}/deepTools/multiBamSummary/{duplicates}/results.npz"
    output:
        "{assayID}/{outdir}/{reference_version}/deepTools/plotCorrelation/{duplicates}/heatmap_SpearmanCorr_readCounts.png",
        "{assayID}/{outdir}/{reference_version}/deepTools/plotCorrelation/{duplicates}/heatmap_SpearmanCorr_readCounts.tab"
    shell:
        """
            {params.deepTools_dir}/plotCorrelation --corData {input.npz} \
                                                   --corMethod spearman \
                                                   --skipZeros \
                                                   --plotTitle "{params.plotTitle}" \
                                                   --whatToPlot heatmap \
                                                   --colorMap RdYlBu \
                                                   --plotNumbers \
                                                   -o {output[0]} \
                                                   --outFileCorMatrix {output[1]}
        """

rule plotPCA:
    params:
        deepTools_dir = home + config["deepTools_dir"],
        plotTitle = lambda wildcards: "PCA - " + wildcards.duplicates
    input:
        npz = "{assayID}/{outdir}/{reference_version}/deepTools/multiBamSummary/{duplicates}/results.npz"
    output:
        "{assayID}/{outdir}/{reference_version}/deepTools/plotPCA/{duplicates}/PCA_readCounts.png"
    shell:
        """
            {params.deepTools_dir}/plotPCA --corData {input.npz} \
                                           --plotFile {output} \
                                           --plotTitle "{params.plotTitle}"
        """

rule bamPEFragmentSize:
    params:
        deepTools_dir = home + config["deepTools_dir"],
        plotTitle = lambda wildcards: "BAM PE " + wildcards.duplicates + " fragment size",
        labels = get_sample_labels
    threads:
        lambda wildcards: int(str(config["program_parameters"]["deepTools"]["threads"]).strip("['']"))
    input:
        expand("{assayID}/{runID}/{outdir}/{reference_version}/bowtie2/duplicates_marked/{unit}.Q{qual}.sorted.MkDup.{suffix}",
               assayID = "ChIP-Seq",
               runID = "SN501_0087_DTremethick_JCSMR_MCF10A_ChIPSeq",
               outdir = config["processed_dir"],
               reference_version = config["references"][REF_GENOME]["version"][0],
               unit = config["samples"]["ChIP-Seq"]["SN501_0087_DTremethick_JCSMR_MCF10A_ChIPSeq"],
               qual = config["alignment_quality"],
               suffix = ["bam"]),
        expand("{assayID}/{runID}/{outdir}/{reference_version}/bowtie2/duplicates_marked/{unit}.Q{qual}.sorted.MkDup.{suffix}",
               assayID = "ChIP-Seq",
               runID = "NB501086_0086_DSTremethick_JCSMR_MCF10A_ChIPseq",
               outdir = config["processed_dir"],
               reference_version = config["references"][REF_GENOME]["version"][0],
               unit = config["samples"]["ChIP-Seq"]["NB501086_0086_DSTremethick_JCSMR_MCF10A_ChIPseq"],
               qual = config["alignment_quality"],
               suffix = ["bam"])
    output:
        "{assayID}/{outdir}/{reference_version}/deepTools/bamPEFragmentSize/{duplicates}/histogram_duplicates_marked.png"
    shell:
        """
            {params.deepTools_dir}/bamPEFragmentSize --bamfiles {input} \
                                                     --samplesLabel {params.labels} \
                                                     --numberOfProcessors {threads} \
                                                     --histogram {output}
        """

rule bamPEFragmentSize_deduplicated:
    params:
        deepTools_dir = home + config["deepTools_dir"],
        plotTitle = lambda wildcards: "BAM PE " + wildcards.duplicates + " fragment size",
        labels = get_sample_labels
    threads:
        lambda wildcards: int(str(config["program_parameters"]["deepTools"]["threads"]).strip("['']"))
    input:
        expand("{assayID}/{runID}/{outdir}/{reference_version}/bowtie2/duplicates_removed/{unit}.Q{qual}.sorted.DeDup.{suffix}",
               assayID = "ChIP-Seq",
               runID = "SN501_0087_DTremethick_JCSMR_MCF10A_ChIPSeq",
               outdir = config["processed_dir"],
               reference_version = config["references"][REF_GENOME]["version"][0],
               unit = config["samples"]["ChIP-Seq"]["SN501_0087_DTremethick_JCSMR_MCF10A_ChIPSeq"],
               qual = config["alignment_quality"],
               suffix = ["bam"]),
        expand("{assayID}/{runID}/{outdir}/{reference_version}/bowtie2/duplicates_removed/{unit}.Q{qual}.sorted.DeDup.{suffix}",
               assayID = "ChIP-Seq",
               runID = "NB501086_0086_DSTremethick_JCSMR_MCF10A_ChIPseq",
               outdir = config["processed_dir"],
               reference_version = config["references"][REF_GENOME]["version"][0],
               unit = config["samples"]["ChIP-Seq"]["NB501086_0086_DSTremethick_JCSMR_MCF10A_ChIPseq"],
               qual = config["alignment_quality"],
               suffix = ["bam"])
    output:
        "{assayID}/{outdir}/{reference_version}/deepTools/bamPEFragmentSize/{duplicates}/histogram_duplicates_removed.png"
    shell:
        """
            {params.deepTools_dir}/bamPEFragmentSize --bamfiles {input} \
                                                     --samplesLabel {params.labels} \
                                                     --numberOfProcessors {threads} \
                                                     --plotTitle "{params.plotTitle}" \
                                                     --histogram {output}
        """

rule plotFingerprint:
    params:
        deepTools_dir = home + config["deepTools_dir"],
        plotTitle = lambda wildcards: "BAM PE " + wildcards.duplicates + " fingerprint",
        labels = get_sample_labels
    threads:
        lambda wildcards: int(str(config["program_parameters"]["deepTools"]["threads"]).strip("['']"))
    input:
        expand("{assayID}/{runID}/{outdir}/{reference_version}/bowtie2/duplicates_marked/{unit}.Q{qual}.sorted.MkDup.{suffix}",
               assayID = "ChIP-Seq",
               runID = "SN501_0087_DTremethick_JCSMR_MCF10A_ChIPSeq",
               outdir = config["processed_dir"],
               reference_version = config["references"][REF_GENOME]["version"][0],
               unit = config["samples"]["ChIP-Seq"]["SN501_0087_DTremethick_JCSMR_MCF10A_ChIPSeq"],
               qual = config["alignment_quality"],
               suffix = ["bam"]),
        expand("{assayID}/{runID}/{outdir}/{reference_version}/bowtie2/duplicates_marked/{unit}.Q{qual}.sorted.MkDup.{suffix}",
               assayID = "ChIP-Seq",
               runID = "NB501086_0086_DSTremethick_JCSMR_MCF10A_ChIPseq",
               outdir = config["processed_dir"],
               reference_version = config["references"][REF_GENOME]["version"][0],
               unit = config["samples"]["ChIP-Seq"]["NB501086_0086_DSTremethick_JCSMR_MCF10A_ChIPseq"],
               qual = config["alignment_quality"],
               suffix = ["bam"])
    output:
        "{assayID}/{outdir}/{reference_version}/deepTools/plotFingerprint/{duplicates}/fingerprints_duplicates_marked.png"
    shell:
        """
            {params.deepTools_dir}/plotFingerprint --bamfiles {input} \
                                                   --numberOfProcessors {threads} \
                                                   --centerReads \
                                                   --plotTitle "{params.plotTitle}" \
                                                   --labels {params.labels} \
                                                   --skipZeros \
                                                   --plotFile {output}
        """

rule plotFingerprint_deduplicated:
    params:
        deepTools_dir = home + config["deepTools_dir"],
        plotTitle = lambda wildcards: "BAM PE " + wildcards.duplicates + " fingerprint",
        labels = get_sample_labels
    threads:
        lambda wildcards: int(str(config["program_parameters"]["deepTools"]["threads"]).strip("['']"))
    input:
        expand("{assayID}/{runID}/{outdir}/{reference_version}/bowtie2/duplicates_removed/{unit}.Q{qual}.sorted.DeDup.{suffix}",
               assayID = "ChIP-Seq",
               runID = "SN501_0087_DTremethick_JCSMR_MCF10A_ChIPSeq",
               outdir = config["processed_dir"],
               reference_version = config["references"][REF_GENOME]["version"][0],
               unit = config["samples"]["ChIP-Seq"]["SN501_0087_DTremethick_JCSMR_MCF10A_ChIPSeq"],
               qual = config["alignment_quality"],
               suffix = ["bam"]),
        expand("{assayID}/{runID}/{outdir}/{reference_version}/bowtie2/duplicates_removed/{unit}.Q{qual}.sorted.DeDup.{suffix}",
               assayID = "ChIP-Seq",
               runID = "NB501086_0086_DSTremethick_JCSMR_MCF10A_ChIPseq",
               outdir = config["processed_dir"],
               reference_version = config["references"][REF_GENOME]["version"][0],
               unit = config["samples"]["ChIP-Seq"]["NB501086_0086_DSTremethick_JCSMR_MCF10A_ChIPseq"],
               qual = config["alignment_quality"],
               suffix = ["bam"])
    output:
        "{assayID}/{outdir}/{reference_version}/deepTools/plotFingerprint/{duplicates}/fingerprints_duplicates_removed.png"
    shell:
        """
            {params.deepTools_dir}/plotFingerprint --bamfiles {input} \
                                                   --numberOfProcessors {threads} \
                                                   --centerReads \
                                                   --plotTitle "{params.plotTitle}" \
                                                   --labels {params.labels} \
                                                   --skipZeros \
                                                   --plotFile {output}
        """
