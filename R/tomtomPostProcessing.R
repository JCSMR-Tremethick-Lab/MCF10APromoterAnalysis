tomtom <- data.table::fread("/home/sebastian/Data/Collaborations/FSU/PromoterSeqCap/SmallFragments/tomtom/cd/TOTALcombined_shH2AZ_Inp_000-125/tomtom.tsv")

summary(tomtom$`q-value`)
summary(tomtom$`E-value`)

length(unique(tomtom[`E-value` < 0.1][`q-value` < 0.1]$Target_ID))
tomtom$extGene <- unlist(lapply(strsplit(tomtom$Target_ID, "_"), function(x) x[1]))
setkey(tomtom, extGene)
rTshH2AZ[tomtom]

rTshH2AZ[which(rTshH2AZ$target_id %in% unique(tomtom[`E-value` < 0.1][`q-value` < 0.1]$extGene))]

vp2  <- ggplot(data = rTshH2AZ[rTshH2AZ[which(rTshH2AZ$target_id %in% tomtom[`E-value` < 0.1]$extGene)]], 
               aes(x = b, y = -(log10(qval)))) +
  xlim(c(-6.5,6.5)) +
  labs(x = "log2 fold-change") +
  geom_point(aes(color = target_id), size = 3) +
  geom_hline(yintercept = 2, linetype = "dotted") +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  theme(legend.position = "none")
vp2

table(rTshH2AZ[which(rTshH2AZ$target_id %in% unique(tomtom[`E-value` < 0.1][`q-value` < 0.1]$extGene))][qval < 0.1]$b > 0)
table(rTshH2AZ[which(rTshH2AZ$target_id %in% unique(tomtom[`E-value` < 0.1][`q-value` < 0.1]$extGene))]$qval < 0.1)


# dot plots of single genes -----------------------------------------------
motif1 <- ggplot(data = kTGenes[tomtom[Query_ID == "YSATTGGC"]$extGene][grep("MCF10AshZD8|MCF10AD8", sample)], 
                     aes(x = condition, y = log2(tpm + 1))) + 
  geom_jitter(aes(color = sample), width = 0.05) + 
  facet_wrap(~target_id) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.4, size = 14))
motif1

motif2 <- ggplot(data = kTGenes[tomtom[Query_ID == "GCCAATCA"]$extGene][grep("MCF10AshZD8|MCF10AD8", sample)], 
                 aes(x = condition, y = log2(tpm + 1))) + 
  geom_jitter(aes(color = sample), width = 0.05) + 
  facet_wrap(~target_id) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.4, size = 14))
motif2

motif3 <- ggplot(data = kTGenes[tomtom[Query_ID == "GGAGGCGG"]$extGene][grep("MCF10AshZD8|MCF10AD8", sample)], 
                 aes(x = condition, y = log2(tpm + 1))) + 
  geom_jitter(aes(color = sample), width = 0.05) + 
  facet_wrap(~target_id) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.4, size = 14))
motif3

motif4 <- ggplot(data = kTGenes[tomtom[Query_ID == "TCCCAGCA"]$extGene][grep("MCF10AshZD8|MCF10AD8", sample)], 
                 aes(x = condition, y = log2(tpm + 1))) + 
  geom_jitter(aes(color = sample), width = 0.05) + 
  facet_wrap(~target_id) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.4, size = 14))
motif4

motif5 <- ggplot(data = kTGenes[tomtom[Query_ID == "TGASTCAB"]$extGene][grep("MCF10AshZD8|MCF10AD8", sample)], 
                 aes(x = condition, y = log2(tpm + 1))) + 
  geom_jitter(aes(color = sample), width = 0.05) + 
  facet_wrap(~target_id) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.4, size = 14))
motif5

# against background of all genes
vp1 <- ggplot(data = rTshH2AZ[qval < 0.1], 
              aes(x = b, y = -(log10(qval)))) +
  geom_point(color = "grey") +
  #  geom_point(aes(color = target_id)) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0)
vp1 + geom_point(data = rTshH2AZ[rTshH2AZ[which(rTshH2AZ$target_id %in% tomtom[`E-value` < 0.1]$extGene)]], aes(color = target_id), size = 4)
