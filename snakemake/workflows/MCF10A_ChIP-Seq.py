_author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2016-12-05"

from snakemake.exceptions import MissingInputException
import os

rule:
    version: 0.1

localrules:
    all, run_kallisto, run_STAR, run_htseq, run_cutadapt

home = os.environ['HOME']

wrapper_dir = home + "/Development/snakemake-wrappers/bio"

include_prefix = home + "/Development/JCSMR-Tremethick-Lab/Breast/snakemake/rules/"

include:
   include_prefix + "perform_cutadapt.py"
include:
    include_prefix + "run_bowtie2.py"
include:
    include_prefix + "bam_processing.py"

rule run_cutadapt:
    input:
        expand("{assayID}/{runID}/{outdir}/{trim_data}/{unit}_{suffix}.QT.CA.fastq.gz",
               assayID = "ChIP-Seq",
               runID = "SN501_0087_DTremethick_JCSMR_MCF10A_ChIPSeq",
               outdir = config["processed_dir"],
               trim_data = config["trim_dir"],
               unit = config["samples"]["ChIP-Seq"]["SN501_0087_DTremethick_JCSMR_MCF10A_ChIPSeq"],
               suffix = ["R1_001", "R2_001"]),
        expand("{assayID}/{runID}/{outdir}/{trim_data}/{unit}_{suffix}.QT.CA.fastq.gz",
               assayID = "ChIP-Seq",
               runID = "NB501086_0086_DSTremethick_JCSMR_MCF10A_ChIPseq",
               outdir = config["processed_dir"],
               trim_data = config["trim_dir"],
               unit = config["samples"]["ChIP-Seq"]["NB501086_0086_DSTremethick_JCSMR_MCF10A_ChIPseq"],
               suffix = ["R1_001", "R2_001"])

rule all:
    input:
        expand("{assayID}/{runID}/{outdir}/{reference_version}/bowtie2/duplicates_marked/{unit}.Q{qual}.sorted.MkDup.{suffix}",
               assayID = "ChIP-Seq",
               runID = "SN501_0087_DTremethick_JCSMR_MCF10A_ChIPSeq",
               outdir = config["processed_dir"],
               reference_version = config["references"]["hg19"]["version"][0],
               unit = config["samples"]["ChIP-Seq"]["SN501_0087_DTremethick_JCSMR_MCF10A_ChIPSeq"],
               qual = config["alignment_quality"],
               suffix = ["bam", "bam.bai"]),
        expand("{assayID}/{runID}/{outdir}/{reference_version}/bowtie2/duplicates_removed/{unit}.Q{qual}.sorted.DeDup.{suffix}",
               assayID = "ChIP-Seq",
               runID = "SN501_0087_DTremethick_JCSMR_MCF10A_ChIPSeq",
               outdir = config["processed_dir"],
               reference_version = config["references"]["hg19"]["version"][0],
               unit = config["samples"]["ChIP-Seq"]["SN501_0087_DTremethick_JCSMR_MCF10A_ChIPSeq"],
               qual = config["alignment_quality"],
               suffix = ["bam", "bam.bai"]),
        expand("{assayID}/{runID}/{outdir}/{reference_version}/bowtie2/duplicates_marked/{unit}.Q{qual}.sorted.MkDup.{suffix}",
               assayID = "ChIP-Seq",
               runID = "NB501086_0086_DSTremethick_JCSMR_MCF10A_ChIPseq",
               outdir = config["processed_dir"],
               reference_version = config["references"]["hg19"]["version"][0],
               unit = config["samples"]["ChIP-Seq"]["NB501086_0086_DSTremethick_JCSMR_MCF10A_ChIPseq"],
               qual = config["alignment_quality"],
               suffix = ["bam", "bam.bai"]),
        expand("{assayID}/{runID}/{outdir}/{reference_version}/bowtie2/duplicates_removed/{unit}.Q{qual}.sorted.DeDup.{suffix}",
               assayID = "ChIP-Seq",
               runID = "NB501086_0086_DSTremethick_JCSMR_MCF10A_ChIPseq",
               outdir = config["processed_dir"],
               reference_version = config["references"]["hg19"]["version"][0],
               unit = config["samples"]["ChIP-Seq"]["NB501086_0086_DSTremethick_JCSMR_MCF10A_ChIPseq"],
               qual = config["alignment_quality"],
               suffix = ["bam", "bam.bai"])
