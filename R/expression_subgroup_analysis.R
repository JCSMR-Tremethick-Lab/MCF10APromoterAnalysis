
# WT vs TGFb --------------------------------------------------------------
design <- model.matrix(~ condition, data = s2c.D6)
so.gene.TGFbD6 <- sleuth::sleuth_prep(s2c.D6, ~ condition,
                               target_mapping = t2gHsap, 
                               aggregation_column = "external_gene_name",
                               filter_fun = filter_function,
                               transformation_function = log2_transform, 
                               max_bootstrap = 30, 
                               read_bootstrap_tpm = T, 
                               extra_bootstrap_summary = T)
so.gene.TGFbD6 <- sleuth::sleuth_fit(so.gene.TGFbD6, formula = design)
so.gene.TGFbD6 <- sleuth::sleuth_fit(so.gene.TGFbD6, ~1, "reduced")
so.gene <- sleuth::sleuth_lrt(so.gene, "reduced", "full")
for (i in colnames(design)[grep("Intercept", colnames(design), invert = T)]){
  so.gene.TGFbD6 <- sleuth::sleuth_wt(so.gene.TGFbD6, i)  
}


# WT vs shZ ---------------------------------------------------------------
s2c.D8 <- s2c[grep("D8", sample)]
s2c.D8$condition <- droplevels(s2c.D8$condition)

design <- model.matrix(~ condition, data = s2c.D8)
so.gene.shZD8 <- sleuth::sleuth_prep(s2c.D8, ~ condition,
                               target_mapping = t2gHsap, 
                               aggregation_column = "external_gene_name",
                               filter_fun = filter_function,
                               transformation_function = log2_transform, 
                               max_bootstrap = 30, 
                               read_bootstrap_tpm = T, 
                               extra_bootstrap_summary = T)
so.gene.shZD8 <- sleuth::sleuth_fit(so.gene.shZD8, formula = design)
so.gene.shZD8 <- sleuth::sleuth_fit(so.gene.shZD8, ~1, "reduced")
so.gene.shZD8 <- sleuth::sleuth_lrt(so.gene.shZD8, "reduced", "full")
for (i in colnames(design)[grep("Intercept", colnames(design), invert = T)]){
  so.gene.shZD8 <- sleuth::sleuth_wt(so.gene.shZD8, i)  
}


# WT vs Ca1a --------------------------------------------------------------
s2c.wt_vs_cancer <- s2c[grep("wt_rep", sample)]
s2c.wt_vs_cancer$condition <- droplevels(s2c.wt_vs_cancer$condition)

design <- model.matrix(~ condition, data = s2c.wt_vs_cancer)
so.gene.wt_vs_cancer <- sleuth::sleuth_prep(s2c.wt_vs_cancer, ~ condition,
                                     target_mapping = t2gHsap, 
                                     aggregation_column = "external_gene_name",
                                     filter_fun = filter_function,
                                     transformation_function = log2_transform, 
                                     max_bootstrap = 30, 
                                     read_bootstrap_tpm = T, 
                                     extra_bootstrap_summary = T)
so.gene.wt_vs_cancer <- sleuth::sleuth_fit(so.gene.wt_vs_cancer, formula = design)
so.gene.wt_vs_cancer <- sleuth::sleuth_fit(so.gene.wt_vs_cancer, ~1, "reduced")
so.gene.wt_vs_cancer <- sleuth::sleuth_lrt(so.gene.wt_vs_cancer, "reduced", "full")
for (i in colnames(design)[grep("Intercept", colnames(design), invert = T)]){
  so.gene.wt_vs_cancer <- sleuth::sleuth_wt(so.gene.wt_vs_cancer, i)  
}


