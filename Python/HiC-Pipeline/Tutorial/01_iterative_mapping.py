import os
import logging
from hiclib import mapping
from mirnylib import h5dict, genome

logging.basicConfig(level=logging.DEBUG)

if not os.path.exists('/home/skurscheid/tmp'):
    os.mkdir('/home/skurscheid/tmp')



# A. Map the reads iteratively.
mapping.iterative_mapping(
    bowtie_path='/usr/local/bin/bowtie2',
    bowtie_index_path='/home/skurscheid/Data/References/Genomes/Homo_sapiens/GRCh38_UCSC/index/hg38',
    fastq_path='SRR027956.sra',
    out_sam_path='SRR027056_1.bam',
    min_seq_len=25,
    len_step=5,
    seq_start=0,
    seq_end=75,
    nthreads=16,  # on intel corei7 CPUs 4 threads are as fast as
                 # 8, but leave some room for you other applications
    #max_reads_per_chunk = 10000000,  #optional, on low-memory machines
    temp_dir='/home/skurscheid/tmp',
    bowtie_flags='--very-sensitive',
    bash_reader='/home/skurscheid/bin/fastq-dump -Z')

mapping.iterative_mapping(
    bowtie_path='/usr/local/bin/bowtie2',
    bowtie_index_path='/home/skurscheid/Data/References/Genomes/Homo_sapiens/GRCh38_UCSC/index/hg38',
    fastq_path='SRR027956.sra',
    out_sam_path='SRR027056_2.bam',
    min_seq_len=25,
    len_step=5,
    seq_start=76,
    seq_end=151,
    nthreads=16,
    #max_reads_per_chunk = 10000000,
    temp_dir='/home/skurscheid/tmp',
    bowtie_flags='--very-sensitive',
    bash_reader='/home/skurscheid/bin/fastq-dump -Z')

# B. Parse the mapped sequences into a Python data structure,
#    assign the ultra-sonic fragments to restriction fragments.
mapped_reads = h5dict.h5dict('mapped_reads.hdf5')
genome_db    = genome.Genome('/home/skurscheid/Data/References/Genomes/Homo_sapiens/GRCh38_UCSC/fasta', readChrms=['#', 'X', 'Y', 'M'])

# continue from here 2016-09-28
mapping.parse_sam(
    sam_basename1='SRR027056_1.bam',
    sam_basename2='SRR027056_2.bam',
    out_dict=mapped_reads,
    genome_db=genome_db,
    enzyme_name='HindIII')
