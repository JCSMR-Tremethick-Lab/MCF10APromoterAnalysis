require(GenomicRanges)
require(clusterProfiler)
require(ComplexHeatmap)

# promoterSeqCap_deepTools_postProcessing.R
mat <- read.table("~/Data/Tremethick/Breast/PromoterSeqCap/processed_data/hg38/deepTools/computeMatrix_referencePoint/EMT_markers.MCF10A_H2AZ_H.matrix.gz", header = F, comment.char = "@")
gr <- GRanges(mat$V1, IRanges(start = mat$V2, end = mat$V3, names = mat$V4), strand = mat$V6)

rownames(mat) <- mat$V4
mat <- mat[,7:ncol(mat)]
mat <- as.matrix(mat)
Heatmap(mat, cluster_columns = F, clustering_method_rows = "centroid", clustering_distance_rows = "pearson")

