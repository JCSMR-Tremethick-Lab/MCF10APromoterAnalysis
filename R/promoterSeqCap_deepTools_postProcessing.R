require(GenomicRanges)
require(clusterProfiler)
require(ComplexHeatmap)
require(fastcluster)
require(data.table)

# promoterSeqCap_deepTools_postProcessing.R
mat <- read.table("~/Data/Tremethick/Breast/PromoterSeqCap/processed_data/hg38/deepTools/computeMatrix_referencePoint/EMT_markers.MCF10A_H2AZ_H.matrix.gz", header = F, comment.char = "@")
gr <- GRanges(mat$V1, IRanges(start = mat$V2, end = mat$V3, names = mat$V4), strand = mat$V6)

rownames(mat) <- mat$V4
mat <- mat[,7:ncol(mat)]
mat <- as.matrix(mat)
Heatmap(mat, cluster_columns = F, clustering_method_rows = "centroid", clustering_distance_rows = "pearson")

mat <- deepToolsUtils::computeMatrixLoader('~/Data/Tremethick/Breast/PromoterSeqCap/processed_data/hg38/deepTools/computeMatrix_referencePoint/seqCapTargets_hg38.MCF10A_H2AZ_L.matrix.gz')

mat$computeMatrix <- mat$computeMatrix[,-1]

mat1 <- as.matrix(mat$computeMatrix)

mat1 <- scale(mat1)
dim(mat1)


