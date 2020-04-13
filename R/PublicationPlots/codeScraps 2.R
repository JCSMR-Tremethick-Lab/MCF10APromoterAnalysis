library(RColorBrewer)
data(genes)
cols <- c(rep("grey80", 24), brewer.pal("YlOrRd", n = 9))
genes$chrom <- factor(genes$chrom, levels=c(paste("chr", 1:22, sep=""), "chrX", "chrY"))
ggparallel(list("path", "chrom"), text.offset=c(0.03, 0,-0.03),
           data = genes,  width=0.1, order=c(1,0), text.angle=0,
           color="white",
           factorlevels =  c(sapply(unique(genes$chrom), as.character),
                             unique(genes$path))) +
  scale_fill_manual(values = cols, guide="none") +
  scale_colour_manual(values = cols, guide="none") +
  coord_flip()

kT1 <- data.table(results$MCF10A_vs_shZ$kallisto_table_genes)
kT1 <- melt(kT1, id.vars = "external_gene_name", measure.vars = c("MCF10AD8_2", "MCF10AD8_3", "MCF10AshZD8_1", "MCF10AshZD8_2", "MCF10AshZD8_3"))
kT1$group <- unlist(lapply(strsplit(as.character(kT1$variable), "_"), function(x) x[1]))
kT1$value <- log2(kT1$value + 1)
setkey(kT1, "external_gene_name")
ggplot(data = kT1[c("H2AFZ", "CBX1", "CBX5")], aes(x = group, y = value)) + geom_jitter(aes(fill = group, colour = group)) + facet_wrap(~external_gene_name) + ylab("log2(cpm + 1)") + theme(axis.text.x = element_text(angle = 45, vjust = 0.7))


ucscTX <- biomaRt::getBM(attributes = c("external_gene_name",
                                        "description",
                                        "entrezgene",
                                        "ucsc"),
                            filters = "ucsc",
                            values = ucscTranscripts$TXNAME,
                            mart = mart)

gather(Fig2Categories, key = "condition", value = "category", c("MCF10A_WT.category", "MCF10A_TGFb.category", "MCF10A_shZ.category", "MCF10CA1A.category"))

