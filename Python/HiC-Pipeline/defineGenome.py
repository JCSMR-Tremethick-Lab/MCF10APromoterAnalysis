from mirnylib import h5dict, genome

"""
Here you need to edit the path to the genome files (see mirnylib.genome file or online help for details)
Other components will automatically

Note that names of fasta files should match exactly to the names of sequences in a fasta file
And that each chromosome has exactly the same sequence
"""
allGenomes = {}

def getGenome(name):
    if name in allGenomes:
        return allGenomes[name]
    if name == "hg19":
        genome_db = genome.Genome("../data/hg19")
    elif name == "hg18":
        genome_db = genome.Genome("../data/hg18")
    elif name == "mm9":
        genome_db = genome.Genome("../data/mm9")
    elif name == "mm10":
        genome_db = genome.Genome("../data/mm10")
    elif name == "hg38":
        genome_db = genome.Genome("/home/skurscheid/Data/References/Genomes/Homo_sapiens/GRCh38_UCSC/fasta")
    else:
        raise ValueError("Genome {0} not defined. Edit defineGenome.py and define it".format(name))
    allGenomes[name] = genome_db
    return genome_db
