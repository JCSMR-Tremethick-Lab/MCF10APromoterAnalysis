t1 <- rbind( table("shZ" = FigS2Data$shZ.group), table(FigS2Data[selected]$shZ.group))
t1
chisq.test(t1, simulate.p.value = T)
t1[1,] <- t1[1,] - t1[2,]
t1
chisq.test(t1, simulate.p.value = T)
t2 <- cbind(t1[,1], rowsum(t1[,-1]))
t2 <- cbind(t1[,1], rowSum(t1[,-1]))
t2 <- cbind(t1[,1], rowSums(t1[,-1]))
t2
fisher.test(t2) # tests if cluster1 changes (dependecy if gene is DE or not) are significant in comparison to promoters to all other clusters together
# repeat until 7 p-values, then BH correction
# odds-ratio is effect size (positive/negative)


rownames(t1) <- c('non_de', 'de')
df <- t1 %>%
  as.data.frame() %>%
  rownames_to_column("condition") %>%
  pivot_longer(-condition, names_to = "cluster") %>%
  pivot_wider(names_from = "condition", values_from = "value")
fit <- glm(
  cbind(non_de, de) ~ 0 + cluster, 
  data = df, 
  family = binomial)

fit %>% tidy(exponentiate = TRUE, conf.int = TRUE)
