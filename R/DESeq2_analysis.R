require(GenomicRanges)
require(pheatmap)
require(RColorBrewer)
require(PoiClaClu)
library(ggplot2)

gr.ensGenes <- GenomicRanges::GRanges(seqnames = ensGenes$chromosome_name,
                                      IRanges(start = ensGenes$start_position,
                                              end = ensGenes$end_position,
                                              names = ensGenes$ensembl_gene_id),
                                      strand = c("+", "-")[match(ensGenes$strand, c(1,-1))],
                                      external_gene_name = ensGenes$external_gene_name,
                                      description = ensGenes$description)


htseq_files <- list.files(paste(pathPrefix,
                                "Data/Tremethick/Breast/RNA-Seq/NB501086_0067_RDomaschenz_JCSMR_RNASeq/processed_data", 
                                runConfig$references$version,
                                "HTSeq/count",
                                sep = "/"),
                          full.names = T)
names(htseq_files) <- sample_id
cts <- deepToolsUtils::makeHTSeqCountMatrix(htseq_files)

cts <- cts[rowSums(cts) > 1, ]
cts <- cts[grep("ENSG", rownames(cts)), ]
cts <- cts[, grep("MCF10A", colnames(cts))]
i1 <- intersect(rownames(cts), names(gr.ensGenes))
length(i1)
s2c.mcf10a <- s2c[grep("MCF10A", s2c$sample),]
s2c.mcf10a$condition <- factor(as.character(s2c.mcf10a$condition))
s2c.mcf10a$condition <- relevel(s2c.mcf10a$condition, "MCF10A_wt")
s2c.mcf10a$sample <- factor(as.character(s2c.mcf10a$sample))

rse <- SummarizedExperiment::SummarizedExperiment(assays = cts,
                                                  rowRanges = gr.ensGenes[i1],
                                                  colData = s2c[grep("MCF10A", s2c$sample),])
colSums(assay(rse))
colData(rse)

dds <- DESeqDataSet(rse, design = ~ condition)

rld <- DESeq2::rlog(dds, blind = F)
dds <- estimateSizeFactors(dds)

par( mfrow = c( 1, 2 ) )
plot(log2(counts(dds, normalized=TRUE)[,c(1,4)] + 1),
     pch=16, cex=0.3)
plot(assay(rld)[,c(1,4)],
     pch=16, cex=0.3)

sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( rld$condition, colnames(rld), sep="-" )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

poisd <- PoissonDistance(t(counts(dds)))
samplePoisDistMatrix <- as.matrix( poisd$dd )
rownames(samplePoisDistMatrix) <- paste( rld$condition, colnames(rld), sep="-" )
colnames(samplePoisDistMatrix) <- NULL
pheatmap(samplePoisDistMatrix,
         clustering_distance_rows=poisd$dd,
         clustering_distance_cols=poisd$dd,
         col=colors)

plotPCA(rld, intgroup = c("condition", "sample"))
data <- plotPCA(rld, intgroup = c("condition", "sample"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
ggplot(data, aes(PC1, PC2, color=sample, shape=condition)) + geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  coord_fixed()


# actual DESeq2 analysis --------------------------------------------------
dds <- DESeq2::DESeq(dds)
res <- results(dds, contrast = c("condition", "MCF10A_shZ", "MCF10A_wt"))
res$external_gene_name <- gr.ensGenes[rownames(res)]$external_gene_name
resLFC1 <- results(dds, lfcThreshold=1, contrast = c("condition", "MCF10A_TGFb", "MCF10A_wt"))
resLFC1 <- resLFC1[order(resLFC1$padj),]
table(resLFC1$padj < 0.1)
topGene <- rownames(resLFC1)[which.min(resLFC1$padj)]
topGenes <- rownames(resLFC1)[1:100]
plotMA(res, ylim = c(-6,6))
with(res[i2, ], {
  points(baseMean, log2FoldChange, col="dodgerblue", cex=2, lwd=2)
  text(baseMean, log2FoldChange, gr.ensGenes[i2]$external_gene_name, pos=2, col="dodgerblue")
})
dev.off()
gr.ensGenes[i2]$external_gene_name
