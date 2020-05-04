library(org.Hs.eg.db)
library(clusterProfiler)
library(data.table)

dataDir <- "./alluvial_plots_sensitivity_input/"

tab1 <- data.table::fread(file = file.path(dataDir, 'diffGenes_WT_2_5_shH2AZ_6.csv'), key = 'target_id')
key(tab1)

