import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from mirnylib import genome
from mirnylib import h5dict
from mirnylib import plotting
from hiclib import binnedData

heatmap = "MCF10A_WT-all-HindIII-40k.byChr"
genome_db = genome.Genome('/home/skurscheid/Data/References/Genomes/Homo_sapiens/GRCh38_UCSC/fasta', readChrms=['#', 'X', 'Y', 'M'])

# Read resolution from the dataset.
raw_heatmap = h5dict.h5dict(heatmap, mode='r')
resolution = int(raw_heatmap['resolution'])
resolution = 40000
# Create a binnedData object, load the data.
BD = binnedData.binnedData(resolution, genome_db)
BD.simpleLoad(heatmap, 'HindIII_MCF10A_WT_all')

# Remove the contacts between loci located within the same bin.
BD.removeDiagonal()

# Remove bins with less than half of a bin sequenced.
BD.removeBySequencedCount(0.5)

# Remove 1% of regions with low coverage.
BD.removePoorRegions(cutoff=1)

# Truncate top 0.05% of interchromosomal counts (possibly, PCR blowouts).
BD.truncTrans(high=0.0005)

# Perform iterative correction.
BD.iterativeCorrectWithoutSS()

# Save the iteratively corrected heatmap.
BD.export('HindIII_MCF10A_WT_all', 'MCF10A_WT-all-HindIII-40k.byChr.hdf5')

# Plot the heatmap directly.
plotting.plot_matrix(np.log(BD.dataDict['HindIII_MCF10A_WT_all']))
plt.savefig("MCF10A_WT-all-HindIII-40k.byChr.pdf")
