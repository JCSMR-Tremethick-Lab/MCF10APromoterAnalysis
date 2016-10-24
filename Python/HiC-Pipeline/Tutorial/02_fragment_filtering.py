from mirnylib import genome
from hiclib import fragmentHiC

# Create a HiCdataset object.
genome_db    = genome.Genome('/home/skurscheid/Data/References/Genomes/Homo_sapiens/GRCh38_UCSC/fasta', readChrms=['#', 'X', 'Y', 'M'])
fragments = fragmentHiC.HiCdataset(
    filename='fragment_dataset.hdf5',
    genome=genome_db,
    maximumMoleculeLength=500,
    mode='w')

# Load the parsed reads into the HiCdataset. The dangling-end filter is applied
# at this stage, with maximumMoleculeLength specified at the initiation of the
# object.
fragments.parseInputData(
    dictLike='mapped_reads.hdf5')

fragments.filterRsiteStart(offset=5)
fragments.filterDuplicates()
fragments.filterLarge()
fragments.filterExtreme(cutH=0.005, cutL=0)
fragments.saveHeatmap('heatmap-res-1M.hdf5', resolution=1000000)
fragments.printMetadata(saveTo="statistics.txt")
