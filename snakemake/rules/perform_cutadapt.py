__author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2016-02-27"

# vim: syntax=python tabstop=4 expandtab
# coding: utf-8


"""
Rules for trimming reads with cutadapt
(http://cutadapt.readthedocs.org/en/latest/guide.html#illumina-truseq)

For usage, include this in your workflow.
"""

rule make_meaningful_filenames:
    params:
        raw_data = config["raw_dir"],
    input:
        read1 = lambda wildcards: wildcards["assayID"] + "/" + wildcards["runID"] + "/fastq/" + config["samples"][wildcards["assayID"]][wildcards["runID"]][wildcards["unit"]][0],
        read2 = lambda wildcards: wildcards["assayID"] + "/" + wildcards["runID"] + "/fastq/" + config["samples"][wildcards["assayID"]][wildcards["runID"]][wildcards["unit"]][1],
    output:
        renamed_read1 = "{assayID}/{runID}/fastq/{unit}_R1_001.fastq.gz",
        renamed_read2 = "{assayID}/{runID}/fastq/{unit}_R2_001.fastq.gz"
    shell:
        """
            mv {input.read1} {output.renamed_read1}; mv {input.read2} {output.renamed_read2}
        """


rule cutadapt_pe:
    params:
        trim_params = config["trim_params"],
        trim_data = config["trim_dir"],
        raw_data = config["raw_dir"],
        cutadapt_dir = config["cutadapt_dir"]
    input:
        read1 = rules.make_meaningful_filenames.output.renamed_read1,
        read2 = rules.make_meaningful_filenames.output.renamed_read2
    output:
        trimmed_read1 = "{assayID}/{runID}/{processed_dir}/{trim_data}/{unit}_R1_001.QT.CA.fastq.gz",
        trimmed_read2 = "{assayID}/{runID}/{processed_dir}/{trim_data}/{unit}_R2_001.QT.CA.fastq.gz"
    shell:
        """
            {params.cutadapt_dir}/cutadapt {params.trim_params} \
                                            -o {output.trimmed_read1} \
                                            -p {output.trimmed_read2} \
                                            {input.read1} \
                                            {input.read2}
        """
