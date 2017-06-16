# export differential UCSC data
rt_MCF10A_shZ <- as.data.table(results[["MCF10A_vs_shZ"]][["sleuth_results"]]$conditionMCF10A_shZ)
rt_MCF10A_shZ <- rt_MCF10A_shZ[!is.na(rt_MCF10A_shZ$pval),]
write_csv(rt_MCF10A_shZ, "~/Data/WorkingData/MCF10A_wt_vs_shZ.csv")
epi_mes_MCF10A_shZ <- rt_MCF10A_shZ[ucscTranscriptsSigEMTCells$ucsc, on = "target_id"]
epi_mes_MCF10A_shZ <- epi_mes_MCF10A_shZ[!is.na(epi_mes_MCF10A_shZ$pval),]
epi_mes_MCF10A_shZ <- merge(epi_mes_MCF10A_shZ, ucscTranscriptsSigEMTCells[,c("ucsc", "epi_mes")], by.x="target_id", by.y="ucsc", all.x = T)
write_csv(epi_mes_MCF10A_shZ, "~/Data/WorkingData/MCF10A_wt_vs_shZ_epi_mes.csv")

rt_MCF10A_TGFb <- as.data.table(results[["MCF10A_vs_TGFb"]][["sleuth_results"]]$conditionMCF10A_TGFb)
rt_MCF10A_TGFb <- rt_MCF10A_TGFb[!is.na(rt_MCF10A_TGFb$pval),]
write_csv(rt_MCF10A_TGFb, "~/Data/WorkingData/MCF10A_wt_vs_TGFb.csv")
epi_mes_MCF10A_TGFb <- rt_MCF10A_TGFb[ucscTranscriptsSigEMTCells$ucsc, on = "target_id"]
epi_mes_MCF10A_TGFb <- epi_mes_MCF10A_TGFb[!is.na(epi_mes_MCF10A_TGFb$pval),]
epi_mes_MCF10A_TGFb <- merge(epi_mes_MCF10A_TGFb, ucscTranscriptsSigEMTCells[,c("ucsc", "epi_mes")], by.x="target_id", by.y="ucsc", all.x = T)
write_csv(epi_mes_MCF10A_TGFb, "~/Data/WorkingData/MCF10A_wt_vs_TGFb_epi_mes.csv")

