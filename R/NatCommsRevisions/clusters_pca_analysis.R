library(rtracklayer)
library(deepToolsUtils)
library(ade4)

MCF10A_WT_Input <- deepToolsUtils::computeMatrixLoader("~/mount/gadi/gdata/PromoterSeqCap/processed_data/GRCh37_hg19_UCSC/deepTools/computeMatrix/MCF10A_WT_Input.Total.gz")

MCF10A_WT_Input.matrix <- MCF10A_WT_Input$computeMatrix[,-1]
MCF10A_WT_Input.rownames <- MCF10A_WT_Input$computeMatrix[,1]

pca1 <- dudi.pca(MCF10A_WT_Input.matrix, scannf = T)

pca1$li

