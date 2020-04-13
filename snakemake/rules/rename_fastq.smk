_author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2015-04-22"

from snakemake.exceptions import MissingInputException
import os

rule:
    version: 0.1


rule make_meaningful_filenames:
    input:
        read1 = lambda wildcards: config[wildcards["unit"]][0],
        read2 = lambda wildcards: config[wildcards["unit"]][1],
    output:
        renamed_read1 = "{unit}_R1_001.fastq.gz",
        renamed_read2 = "{unit}_R2_001.fastq.gz"
    shell:
        """
            mv {input.read1} {output.renamed_read1}; mv {input.read2} {output.renamed_read2}
        """

rule rename_fastq:
   input:
       expand("{unit}_{suffix}.fastq.gz",
              unit = ["MCF10A_wt_rep1",
                      "MCF10A_wt_rep2",
                      "MCF10A_shZ_rep1",
                      "MCF10A_shZ_rep2",
                      "MCF10A_TGFb_rep1",
                      "MCF10A_TGFb_rep2",
                      "MCF10Ca1a_wt_rep1",
                      "MCF10Ca1a_wt_rep2",
                      "MCF10Ca1a_shZ_rep1",
                      "MCF10Ca1a_shZ_rep2"],
              suffix = ["R1_001", "R2_001"])


