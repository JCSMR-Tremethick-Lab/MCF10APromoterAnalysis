library(org.Hs.eg.db)
library(clusterProfiler)
library(data.table)
library(msigdbr)
library(bitr)

dataDir <- "./alluvial_plots_sensitivity_input/"
# enricher
m_t2g <- msigdbr(species = "Homo sapiens") %>% dplyr::select(gs_name, human_gene_symbol)
m_t2g_H <- msigdbr(species = "Homo sapiens", category = 'H') %>% dplyr::select(gs_name, human_gene_symbol)
m_t2g_C2 <- msigdbr(species = "Homo sapiens", category = 'C2') %>% dplyr::select(gs_name, human_gene_symbol)
m_t2g_C5 <- msigdbr(species = "Homo sapiens", category = 'C5') %>% dplyr::select(gs_name, human_gene_symbol)


# first enrichment
tab1 <- data.table::fread(file = file.path(dataDir, 'diffGenes_WT_2_5_shH2AZ_6.csv'), key = 'target_id')
key(tab1)

gois <- tab1$target_id
em_tab1 <- enricher(gois, universe = universe, pvalueCutoff = 0.01, TERM2GENE = m_t2g)
em_H_tab1 <- enricher(gois, universe = universe, pvalueCutoff = 0.01, TERM2GENE = m_t2g_H)
em_C2_tab1 <- enricher(gois, universe = universe, pvalueCutoff = 0.01, TERM2GENE = m_t2g_C2)
em_C5_tab1 <- enricher(gois, universe = universe, pvalueCutoff = 0.01, TERM2GENE = m_t2g_C5)

em_C2_tab1@result %>% as.data.table
png(file = file.path(dataDir, "cluster_profiler_sensitivity_input_2_5_wt_shz_6.png"), width = 1024, height = 1024)
dotplot(em_C2_tab1)
dev.off()

# second enrichment
tab2 <- data.table::fread(file = file.path(dataDir, 'diffGenes_WT_2_5_shH2AZ_not_6.csv'), key = 'target_id')
gois <- tab2$target_id

em_tab2 <- enricher(gois, universe = universe, pvalueCutoff = 0.01, TERM2GENE = m_t2g)
em_H_tab2 <- enricher(gois, universe = universe, TERM2GENE = m_t2g_H, qvalueCutoff = 0.9, pvalueCutoff = 0.9)
em_C2_tab2 <- enricher(gois, universe = universe, TERM2GENE = m_t2g_C2, qvalueCutoff = 0.1)
em_C5_tab2 <- enricher(gois, universe = universe, TERM2GENE = m_t2g_C5, qvalueCutoff = 0.9, pvalueCutoff = 0.9)

em_tab2@result %>% as.data.table -> em_tab2_result
em_tab2_result[(em_tab2_result$qvalue < 0.000001)]

png(file = file.path(dataDir, "cluster_profiler_sensitivity_input_2_5_wt_shz_not_6.png"), width = 1024, height = 1024)
dotplot(em_C2_tab2)
dev.off()
