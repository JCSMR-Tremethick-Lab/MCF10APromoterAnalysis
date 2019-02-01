_author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2015-04-22"

from snakemake.exceptions import MissingInputException
import os

rule:
    version: 0.3

localrules:
    all, run_kallisto, run_STAR, run_htseq, run_cutadapt

wrapper_dir = os.environ['HOME'] + "/Development/snakemake-wrappers/bio"

include_prefix= os.environ['HOME'] + "/Development/JCSMR-Tremethick-Lab/Breast/snakemake/rules/"

include:
     include_prefix + "perform_fastqc.py"
include:
    include_prefix + "perform_cutadapt.py"
include:
    include_prefix + "run_kallisto.py"
include:
    include_prefix + "run_STAR.py"

rule run_kallisto:
    input:
        expand("{assayID}/NB501086_0067_RDomaschenz_JCSMR_RNASeq/{outdir}/{reference_version}/kallisto/{unit}",
               assayID = "RNA-Seq",
               outdir = config["processed_dir"],
               reference_version = config["references"]["hg19"]["version"],
               unit = config["samples"]["RNA-Seq"]),
        expand("{assayID}/NB501086_0082_RDomaschenz_JCSMR_mRNAseq/{outdir}/{reference_version}/kallisto/{unit}",
               assayID = "RNA-Seq_run2",
               outdir = config["processed_dir"],
               reference_version = config["references"]["hg19"]["version"],
               unit = config["samples"]["RNA-Seq_run2"])

rule run_STAR:
    input:
        expand("{assayID}/NB501086_0067_RDomaschenz_JCSMR_RNASeq/{outdir}/{reference_version}/STAR/full/{unit}.aligned.bam",
               assayID = "RNA-Seq",
               outdir = config["processed_dir"],
               reference_version = config["references"]["hg19"]["version"],
               unit = config["samples"]["RNA-Seq"]),
        expand("{assayID}/NB501086_0082_RDomaschenz_JCSMR_mRNAseq/{outdir}/{reference_version}/STAR/full/{unit}.aligned.bam",
               assayID = "RNA-Seq_run2",
               outdir = config["processed_dir"],
               reference_version = config["references"]["hg19"]["version"],
               unit = config["samples"]["RNA-Seq_run2"])

rule run_STAR_untrimmed:
    input:
        expand("{assayID}/{runID}/{outdir}/{reference_version}/STAR/full/untrimmed/{unit}.aligned.bam",
               assayID = "RNA-Seq",
               runID = "NB501086_0067_RDomaschenz_JCSMR_RNASeq",
               outdir = config["processed_dir"],
               reference_version = config["references"]["hg19"]["version"],
               unit = config["samples"]["RNA-Seq"]),
        expand("{assayID}/{runID}/{outdir}/{reference_version}/STAR/full/untrimmed/{unit}.aligned.bam",
               assayID = "RNA-Seq_run2",
               runID = "NB501086_0082_RDomaschenz_JCSMR_mRNAseq",
               outdir = config["processed_dir"],
               reference_version = config["references"]["hg19"]["version"],
               unit = config["samples"]["RNA-Seq_run2"])

rule run_htseq:
    input:
        expand("{assayID}/NB501086_0067_RDomaschenz_JCSMR_RNASeq/{outdir}/{reference_version}/HTSeq/count/{unit}.txt",
               assayID = "RNA-Seq",
               outdir = config["processed_dir"],
               reference_version = config["references"]["hg19"]["version"],
               unit = config["samples"]["RNA-Seq"]),
        expand("{assayID}/NB501086_0082_RDomaschenz_JCSMR_mRNAseq/{outdir}/{reference_version}/HTSeq/count/{unit}.txt",
               assayID = "RNA-Seq_run2",
               outdir = config["processed_dir"],
               reference_version = config["references"]["hg19"]["version"],
               unit = config["samples"]["RNA-Seq_run2"])

rule run_cutadapt:
    input:
        expand("{assayID}/{runID}/{outdir}/{trim_data}/{unit}_{suffix}.QT.CA.fastq.gz",
               assayID = "RNA-Seq",
               runID = "NB501086_0067_RDomaschenz_JCSMR_RNASeq",
               outdir = config["processed_dir"],
               trim_data = config["trim_dir"],
               unit = config["samples"]["RNA-Seq"],
               suffix = ["R1_001", "R2_001"]),
        expand("{assayID}/{runID}/{outdir}/{trim_data}/{unit}_{suffix}.QT.CA.fastq.gz",
               assayID = "RNA-Seq_run2",
               runID = "NB501086_0082_RDomaschenz_JCSMR_mRNAseq",
               outdir = config["processed_dir"],
               trim_data = config["trim_dir"],
               unit = config["samples"]["RNA-Seq_run2"],
               suffix = ["R1_001", "R2_001"])

rule rename_fastq:
       input:
           expand("{assayID}/{runID}/fastq/{unit}_{suffix}.fastq.gz",
              assayID = "RNA-Seq",
              runID = "NB501086_0067_RDomaschenz_JCSMR_RNASeq",
              unit = config["samples"]["RNA-Seq"],
              suffix = ["R1_001", "R2_001"]))

rule all:
    input:
        expand("{assayID}/{runID}/{outdir}/{trim_data}/{unit}_{suffix}.QT.CA.fastq.gz",
               assayID = "RNA-Seq",
               runID = "NB501086_0067_RDomaschenz_JCSMR_RNASeq",
               outdir = config["processed_dir"],
               trim_data = config["trim_dir"],
               unit = config["samples"]["RNA-Seq"],
               suffix = ["R1_001", "R2_001"]),
        expand("{assayID}/NB501086_0067_RDomaschenz_JCSMR_RNASeq/{outdir}/{reference_version}/kallisto/{unit}",
               assayID = "RNA-Seq",
               outdir = config["processed_dir"],
               reference_version = config["references"]["hg19"]["version"],
               unit = config["samples"]["RNA-Seq"]),
        expand("{assayID}/NB501086_0067_RDomaschenz_JCSMR_RNASeq/{outdir}/{reference_version}/STAR/full/{unit}.aligned.bam",
               assayID = "RNA-Seq",
               outdir = config["processed_dir"],
               reference_version = config["references"]["hg19"]["version"],
               unit = config["samples"]["RNA-Seq"]),
        expand("{assayID}/NB501086_0067_RDomaschenz_JCSMR_RNASeq/{outdir}/{reference_version}/HTSeq/count/{unit}.txt",
               assayID = "RNA-Seq",
               outdir = config["processed_dir"],
               reference_version = config["references"]["hg19"]["version"],
               unit = config["samples"]["RNA-Seq"]),
        # second run
        expand("{assayID}/{runID}/{outdir}/{trim_data}/{unit}_{suffix}.QT.CA.fastq.gz",
               assayID = "RNA-Seq_run2",
               runID = "NB501086_0082_RDomaschenz_JCSMR_mRNAseq",
               outdir = config["processed_dir"],
               trim_data = config["trim_dir"],
               unit = config["samples"]["RNA-Seq_run2"],
               suffix = ["R1_001", "R2_001"]),
        expand("{assayID}/NB501086_0082_RDomaschenz_JCSMR_mRNAseq/{outdir}/{reference_version}/kallisto/{unit}",
               assayID = "RNA-Seq_run2",
               outdir = config["processed_dir"],
               reference_version = config["references"]["hg19"]["version"],
               unit = config["samples"]["RNA-Seq_run2"]),
        expand("{assayID}/NB501086_0082_RDomaschenz_JCSMR_mRNAseq/{outdir}/{reference_version}/STAR/full/{unit}.aligned.bam",
               assayID = "RNA-Seq_run2",
               outdir = config["processed_dir"],
               reference_version = config["references"]["hg19"]["version"],
               unit = config["samples"]["RNA-Seq_run2"]),
        expand("{assayID}/NB501086_0082_RDomaschenz_JCSMR_mRNAseq/{outdir}/{reference_version}/HTSeq/count/{unit}.txt",
               assayID = "RNA-Seq_run2",
               outdir = config["processed_dir"],
               reference_version = config["references"]["hg19"]["version"],
               unit = config["samples"]["RNA-Seq_run2"])
