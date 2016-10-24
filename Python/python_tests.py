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

with open("snakemake/configs/config_PromoterCapSeq.json") as data_file:
    config = json.load(data_file)
