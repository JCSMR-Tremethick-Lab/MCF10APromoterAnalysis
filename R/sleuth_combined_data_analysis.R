s2c$condition <- factor(unlist(lapply(strsplit(as.character(condition), "D"), function(x) x[1])))
s2c$timepoint <- factor(paste("D", unlist(lapply(strsplit(unlist(lapply(strsplit(as.character(condition), "D"), function(x) x[2])), "_"), function(x) x[1])), sep = ""))

design <- model.matrix(~ condition, data = s2c)
print(paste("Processing ", x, " transcript-level analysis",sep = ""))
#-----------------------------------------------------------------------------
# transcript-level DE
so <- sleuth::sleuth_prep(s2c, ~ condition, target_mapping = t2g, max_bootstrap = 30, extra_bootstrap_summary = T)
so <- sleuth::sleuth_fit(so, formula = design)
so <- sleuth::sleuth_fit(so, ~1, "reduced")
so <- sleuth::sleuth_lrt(so, "reduced", "full")
for (i in colnames(design)[grep("Intercept", colnames(design), invert = T)]){
  so <- sleuth::sleuth_wt(so, i)  
}
rt.list <- lapply(colnames(design)[grep("Intercept", colnames(design), invert = T)], function(x){
  rt <- sleuth::sleuth_results(so, x)
  rt <- rt[order(rt$qval),]
})
names(rt.list) <- colnames(design)[grep("Intercept", colnames(design), invert = T)]
kt <- sleuth::kallisto_table(so, normalized = T, include_covariates = T)
kt_wide <- tidyr::spread(kt[, c("target_id", "sample", "tpm")], sample, tpm)
kt_wide <- data.table::as.data.table(kt_wide)
data.table::setkey(kt_wide, "target_id")
so.gene <- sleuth::sleuth_prep(s2c, ~ condition, target_mapping = t2g, aggregation_column = "ext_gene", max_bootstrap = 30, extra_bootstrap_summary = T)
so.gene <- sleuth::sleuth_fit(so.gene, formula = design)
so.gene <- sleuth::sleuth_fit(so.gene, ~1, "reduced")
so.gene <- sleuth::sleuth_lrt(so.gene, "reduced", "full")

# combined
design.combined <- model.matrix(~ condition + timepoint, data = s2c)
so.combined <- sleuth::sleuth_prep(s2c, design.combined, target_mapping = t2g, max_bootstrap = 30, extra_bootstrap_summary = T)
so.combined <- sleuth::sleuth_fit(so.combined, formula = design.combined)
so.combined <- sleuth::sleuth_fit(so.combined, ~1, "reduced")
so.combined <- sleuth::sleuth_lrt(so.combined, "reduced", "full")
for (i in colnames(design.combined)[grep("condition", colnames(design.combined), invert = F)]){
  so.combined <- sleuth::sleuth_wt(so.combined, i)  
}

