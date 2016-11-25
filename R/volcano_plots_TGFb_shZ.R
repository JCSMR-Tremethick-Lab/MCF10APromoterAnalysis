# volcano plot - TGFb vs WT
emtGenes <- resultsCompressed[[2]]$sleuth_results_genes[["conditionMCF10A_TGFb"]][which(resultsCompressed[[2]]$sleuth_results_genes[["conditionMCF10A_TGFb"]]$target_id %in% hsap.qPCRGenesTab$ensembl_gene_id),]
emtGenes$FC_estimated <- log2(exp(emtGenes$b))
emtGenes <- merge(emtGenes, ensGenes[,c("ensembl_gene_id", "external_gene_name")], by.x = "target_id", by.y = "ensembl_gene_id")
pdf("Volcano_plot_EMT_genes_TGFb.pdf")
plot(emtGenes$b, -log10(emtGenes$qval), axes = F, xlab = "", ylab = "", frame = F, xlim = c(-10,10), cex = 0.3, pch = 16)
points(emtGenes[which(-log10(emtGenes$qval) >= 1), "b"], -log10(emtGenes[which(-log10(emtGenes$qval) >= 1), "qval"]), col = "red", pch = 16, cex = 1.2)
axis(2, pos = 0, lwd = 3) #, at = c(seq(0,200,2)), labels = c("", seq(2,200,2)))
axis(1, pos = 0, lwd = 3)
mtext("-log10(q-value)", side = 2)
mtext("beta value", side = 1, line = 2)
abline(h = 1, col = "red", lty = 2, lwd = 2)
# only adding labels to genes with adjuste p-value <= 0.01
text(emtGenes[which(-log10(emtGenes$qval) >= 2), "b"], -log10(emtGenes[which(-log10(emtGenes$qval) >= 2), "qval"]), 
     labels = emtGenes[which(-log10(emtGenes$qval) >= 2), "external_gene_name" ],
     cex = 0.7,
     pos = 4, offset = 0.3)
text(-4,1.2, "< 0.1 [adjusted p-value]", cex = 0.7)
dev.off()

# volcano plot - shZ vs WT
emtGenes.shZ <- resultsCompressed[[1]]$sleuth_results_genes[["conditionMCF10A_shZ"]][which(resultsCompressed[[1]]$sleuth_results_genes[["conditionMCF10A_shZ"]]$target_id %in% hsap.qPCRGenesTab$ensembl_gene_id),]
emtGenes.shZ$FC_estimated <- log2(exp(emtGenes.shZ$b))
emtGenes.shZ <- merge(emtGenes.shZ, ensGenes[,c("ensembl_gene_id", "external_gene_name")], by.x = "target_id", by.y = "ensembl_gene_id")
pdf("Volcano_plot_EMT_genes_shZ.pdf")
plot(emtGenes.shZ$b, -log10(emtGenes.shZ$qval), axes = F, xlab = "", ylab = "", frame = F, xlim = c(-10,10), cex = 0.3, pch = 16)
points(emtGenes.shZ[which(-log10(emtGenes.shZ$qval) >= 1), "b"], -log10(emtGenes.shZ[which(-log10(emtGenes.shZ$qval) >= 1), "qval"]), col = "red", pch = 16, cex = 1.2)
axis(2, pos = 0, lwd = 3) #, at = c(seq(0,20,2)), labels = c("", seq(2,20,2)))
axis(1, pos = 0, lwd = 3)
mtext("-log10(q-value)", side = 2)
mtext("beta value", side = 1, line = 2)
abline(h = 1, col = "red", lty = 2, lwd = 2)
# only adding labels to genes with adjuste p-value <= 0.01
text(emtGenes.shZ[which(-log10(emtGenes.shZ$qval) >= 2), "b"], -log10(emtGenes.shZ[which(-log10(emtGenes.shZ$qval) >= 2), "qval"]), 
     labels = emtGenes.shZ[which(-log10(emtGenes.shZ$qval) >= 2), "external_gene_name" ],
     cex = 0.7,
     pos = 4, offset = 0.3)
text(-4,1.2, "< 0.1 [adjusted p-value]", cex = 0.7)
dev.off()

