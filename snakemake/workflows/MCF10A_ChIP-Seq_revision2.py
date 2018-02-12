_author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2016-12-05"

from snakemake.exceptions import MissingInputException
import os

rule:
    version: 0.2

localrules:
    all

home = os.environ['HOME']

wrapper_dir = home + "/Development/snakemake-wrappers/bio"

include_prefix = home + "/Development/JCSMR-Tremethick-Lab/Breast/snakemake/rules/"

include:
    include_prefix + "run_bowtie2.py"

# define global variables such as reference version of genome so that it can be accessed
# throughout the whole worfklow
REF_GENOME = config["references"]["genomes"][0]

rule all:
    input:
        expand("{assayID}/{runID}/{outdir}/{reference_version}/bowtie2/{duplicates}/{sample}.Q{qual}.sorted.{suffix}",
               assayID = "ChIP-Seq",
               runID = "SN501_0087_DTremethick_JCSMR_MCF10A_ChIPSeq",
               outdir = config["processed_dir"],
               reference_version = config["references"][REF_GENOME]["version"],
               duplicates = ["duplicates_marked", "duplicates_removed"],
               sample = config["samples"]["ChIP-Seq"]["SN501_0087_DTremethick_JCSMR_MCF10A_ChIPSeq"],
               qual = config["alignment_quality"],
               suffix = ["bam", "bam.bai"]),
        expand("{assayID}/{runID}/{outdir}/{reference_version}/bowtie2/{duplicates}/{sample}.Q{qual}.sorted.{suffix}",
               assayID = "ChIP-Seq",
               runID = "NB501086_0086_DSTremethick_JCSMR_MCF10A_ChIPseq",
               outdir = config["processed_dir"],
               reference_version = config["references"][REF_GENOME]["version"],
               duplicates = ["duplicates_marked", "duplicates_removed"],
               sample = config["samples"]["ChIP-Seq"]["NB501086_0086_DSTremethick_JCSMR_MCF10A_ChIPseq"],
               qual = config["alignment_quality"],
               suffix = ["bam", "bam.bai"]),
