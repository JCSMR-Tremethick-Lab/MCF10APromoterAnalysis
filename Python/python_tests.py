import json
from pprint import pprint

def get_replicates_input(wildcards):
    fn = []
    for i in config["samples"][wildcards["digest"]][wildcards["Input"]][wildcards["sample"] + "_" + wildcards["Input"] + "_" + wildcards["digest"]]:
        fn.append("processed_data/hg38/duplicates_removed/" + i + ".DeDup.sorted.fastq_q20.bam")
    return(fn)

def get_replicates_chip(wildcards):
    fn = []
    for i in config["samples"][wildcards["digest"]][wildcards["ChIP"]][wildcards["sample"] + "_" + wildcards["ChIP"] + "_" + wildcards["digest"]]:
        fn.append("processed_data/hg38/duplicates_removed/" + i + ".DeDup.sorted.fastq_q20.bam")
    return(fn)

wildcards = dict()
wildcards = {"digest" : "H", "sample" : "MCF10A_TGFb", "Input" : "Input", "ChIP" : "H2AZ"}

with open("snakemake/configs/config_ChIP-Seq.json") as data_file:
    config = json.load(data_file)

with open("configs/config.json") as data_file:
    config = json.load(data_file)

with open("/Users/u1001407/Development/JCSMR-Tremethick-Lab/H2AZ_EMT/snakemake/configs/config.json") as data_file:
    config = json.load(data_file)


with open("config.json") as data_file:
    config = json.load(data_file)


' '.join("{!s}={!r}".format(key, val) for (key, val) in cli_parameters(wildcards).items())

# this returns a single string from multiple dictionary items, including keys and values - great for generating program CLI parameters
' '.join("{!s}={!s}".format(key, val) for (key, val) in config["program_parameters"]["deepTools"]["computeMatrix"]["reference-point"].items())

def cli_parameters_computeMatrix(wildcards):
    a = config["program_parameters"][wildcards.application][wildcards.tool][wildcards.mode]
    if wildcards.mode == "reference-point":
        a["--referencePoint"] = wildcards.referencePoint
    return(a)

def cli_parameters_bamCoverage(wildcards):
    a = config["program_parameters"][wildcards["application"]][wildcards["tool"]][wildcards["mode"]]
    b = str()
    for (key, val) in a.items():
        if val == " ":
            f = key + " "
            b = b + f
        else:
            f = key + "=" + val + " "
            b = b + f
    if wildcards["mode"] == "MNase":
        b = b + "--MNase"
    return(b.rstrip())


b = b + join("{!s}={!s}".format(key, val))
    #b = ' '.join("{!s}={!s}".format(key, val.strip("\\'")) for (key, val) in a.items())


cli_parameters(wildcards)
cli_parameters_bamCoverage(wildcards)

', '.join("{!s}={!r}".format(key, val) for (key, val) in config["samples"]["ChIP-Seq"]["NB501086_0011_MNekrasov_MDCK_JCSMR_ChIPseq"].items())
' '.join("{!s}".format(key) for (key) in config["samples"]["ChIP-Seq"]["NB501086_0011_MNekrasov_MDCK_JCSMR_ChIPseq"].keys())

' '.join("{!s}".format(key) for (key) in config["samples"]["ChIP-Seq"].keys())
config["samples"]["ChIP-Seq"][runIDs]

pprint(key,val) for (key, val) in config["samples"]["ChIP-Seq"]["NB501086_0011_MNekrasov_MDCK_JCSMR_ChIPseq"].items())

wildcards = dict()
wildcards = {"mode" : "MNase", "tool" : "bamCoverage", "application" : "deepTools"}
wildcards = {"mode" : "normal", "tool" : "bamCoverage", "application" : "deepTools"}

wildcards = {"referencePoint" : "TSS", "mode" : "reference-point"}
wildcards = {"ChIP-Seq" : "NB501086_0086_DSTremethick_JCSMR_MCF10A_ChIPseq"}
def get_sample_labels(wildcards):
    fn = []
    runIDs = config["samples"]["ChIP-Seq"].keys()
    for i in runIDs:
        for k in config["samples"][wildcards["ChIP-Seq"]][i].keys():
            fn.append(k)
    return(fn)


def bam_merge_input(wildcards):
    fn = []
    path = "/".join((wildcards["assayID"],
                     wildcards["runID"],
                     wildcards["outdir"],
                     wildcards["reference_version"],
                     wildcards["application"],
                     wildcards["duplicates"]))
    for i in config["samples"]["ChIP-Seq"]["replicates"][wildcards["sample"]]:
        fn.append("/".join((path, ".".join((i, "Q20.sorted.bam")))))
    return(fn)

wildcards = {"assayID" : "ChIP-Seq",
             "runID" : "merged",
             "outdir" : "processed_data",
             "reference_version": config["references"][REF_GENOME]["version"][0],
             "application" : "bowtie2",
             "duplicates" : "duplicates_marked",
             "sample" : "H2AZ_10A_high"}

bam_merge_input(wildcards)


expand("{assayID}/{runID}/{outdir}/{reference_version}/bowtie2/{duplicates}/{sample}.Q{qual}.sorted.{suffix}",
   assayID = "ChIP-Seq",
   runID = "SN501_0087_DTremethick_JCSMR_MCF10A_ChIPSeq",
   outdir = config["processed_dir"],
   reference_version = config["references"][REF_GENOME]["version"],
   duplicates = ["duplicates_marked", "duplicates_removed"],
   sample = config["samples"]["ChIP-Seq"]["SN501_0087_DTremethick_JCSMR_MCF10A_ChIPSeq"],
   qual = config["alignment_quality"],
   suffix = ["bam", "bam.bai"])
