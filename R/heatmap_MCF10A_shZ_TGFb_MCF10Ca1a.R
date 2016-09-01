# heatmap of samples using MCF10A_wt as reference -------------------------
df1 <- resultsCompressed[["MCF10A"]]$kallisto_table_genes[,c(1:9)]
df1 <- df1[, c("ensembl_gene_id", s2c.list[["MCF10A"]]$sample)]
df1 <- df1[complete.cases(df1),]
df2 <- log2(df1[, c(2:9)] + 1)
filter <- apply(df2, 1, function(y) length(y[y>2])>=0.1)
df2 <- as.matrix(df2[filter, ])
colnames(df2) <- as.character(s2c.list[["MCF10A"]]$condition)
sd1 <- apply(df2, 1, sd)
length(sd1)
pdf("Heatmap_MCF10A_MCF10AshZ_MCF10ATGFb_MCF10Ca1a.pdf")
heatmap.3(df2[sd1 > 2,], 
          trace = "none",
          cexCol = 0.6,
          labRow = NULL)
dev.off()
