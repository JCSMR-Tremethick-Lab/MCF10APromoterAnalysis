library(sleuth)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(AnnotationDbi)
library(bitr)
library(data.table)
library(tximport)

dataDir <- "./alluvial_plots_sensitivity_h2az/"

# prepare annotation table
t2g <- AnnotationDbi::select(TxDb.Hsapiens.UCSC.hg19.knownGene, 
       keys = keys(TxDb.Hsapiens.UCSC.hg19.knownGene, keytype = 'TXNAME'),
       keytype = 'TXNAME',
       columns = c('TXNAME', 'GENEID'))
setDT(t2g)
t2g <- t2g[!is.na(t2g$GENEID)]
t2g <- merge(t2g, EG2Symbol(unique(t2g$GENEID), organism = 'human'), by.x = 'GENEID', by.y = 'entrezgene', all.x = T, all.y = F)
setDT(t2g)
setnames(t2g, 'TXNAME', 'target_id')

# load kallisto data
kallistoDir <- '~/Data/Tremethick/Breast/RNA-Seq/NB501086_0067_RDomaschenz_JCSMR_RNASeq/processed_data/GRCh37_hg19_UCSC_ERCC/kallisto'
samples <- list.files(kallistoDir, pattern = 'wt')

s2c <- data.table(sample = samples, condition = c(rep('WT', 2), rep('ca1a', 2)))
kal_dirs <- sapply(samples, function(id) file.path(kallistoDir, id))
s2c <- dplyr::mutate(s2c, path = kal_dirs)
s2c$condition <- as.factor(s2c$condition)
s2c$condition <- relevel(s2c$condition, ref = 'WT')

design <- model.matrix(~ condition, data = s2c)
files <- file.path(s2c$path, 'abundance.h5')
txi <- tximport(files, type = 'kallisto', txOut = T)

so <- sleuth_prep(s2c, ~ condition, 
                  norm_fun_counts = norm_factors_ercc,
                  norm_fun_tpm = norm_factors_ercc,
                  target_mapping = t2g[,c('target_id', 'hgnc_symbol')], 
                  aggregation_column = 'hgnc_symbol', 
                  extra_bootstrap_summary = F,
                  read_bootstrap_tpm = F,
                  max_bootstrap = 30,
                  gene_mode = T)

so <- sleuth::sleuth_fit(so, formula = design)
so <- sleuth::sleuth_fit(so, ~1, "reduced")
so <- sleuth::sleuth_lrt(so, "reduced", "full")
so <- sleuth::sleuth_wt(so, 'conditionca1a') 

rt <- sleuth::sleuth_results(so, 'conditionca1a', show_all = F)
setDT(rt)
rt[order(b,decreasing=TRUE),]
save(so, file = file.path(dataDir, 'so_wt_ca1a.rda'))

