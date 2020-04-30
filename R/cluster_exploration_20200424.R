library(data.table)
library(deepToolsUtils)
library(magrittr)
library(reshape2)
library(cluster)
library(fpc)
library(TSrepr)
library(clusterCrit)
library(ggplot2)
library(parallel)

source_dir <- '~/mount/gadi/gdata/PromoterSeqCap/processed_data/GRCh37_hg19_UCSC/deepTools/computeMatrix/'
dest_dir <- '~/Data/Collaborations/FSU/PromoterSeqCap/cluster_exploration'

list.files(source_dir, pattern = 'gz')
mcf10_wt_input_total <- deepToolsUtils::computeMatrixLoader(file.path(source_dir, 'MCF10A_WT_Input.Total.gz'))
matrix <- mcf10_wt_input_total$computeMatrix
rownames(matrix) <- mcf10_wt_input_total$computeMatrixRows

matrix_seasprof <- TSrepr::repr_matrix(matrix, func = repr_seas_profile, args = list(freq = 40, func = mean),
                                      normalise = TRUE, func_norm = norm_z)
dim(matrix_seasprof)  

  
matrixNorm <- scale(matrix, center = F, scale = T)
matrixNormPCA <- prcomp(matrix, center = T, scale. = TRUE, rank. = 500)

Fig1Asorts <- data.table::fread('~/Data/Collaborations/FSU/PromoterSeqCap/PublicationFigures/NatComms_revisions/Fig1A_Total_10A_Input_k7sorting.tsv')
l1 <- list()
l1[['Fig1A']] <- Fig1Asorts
l1 <- lapply(l1, function(x) {
  ucscID <- unlist(lapply(strsplit(x$gene, ";"), function(x) x[4]))
  extGene <- unlist(lapply(strsplit(x$gene, ";"), function(x) x[6]))
  x$ucscID <- ucscID
  x$extGene <- extGene
  x$chr <-  unlist(lapply(strsplit(x$gene, ";"), function(x) x[1]))
  x$start <-  unlist(lapply(strsplit(x$gene, ";"), function(x) x[2]))
  x$end <-  unlist(lapply(strsplit(x$gene, ";"), function(x) x[3]))
  return(x)
})
l1
setkey(l1$Fig1A, ucscID)

#original data
s1 <- sample(mcf10_wt_input_total$computeMatrixRows, 5000, replace = F)

orig_data <- data.table(melt(data.table(matrixNorm)))
orig_data[, Window := rep(1:ncol(matrix), each = nrow(matrixNorm))]
orig_data[, Promoter := rep(mcf10_wt_input_total$computeMatrixRows, ncol(matrixNorm))]
setkey(orig_data, Promoter)

orig_data <- merge(orig_data, l1$Fig1A[,c('ucscID', 'group1')], by.x = 'Promoter', by.y = 'ucscID', all.x = T)
colnames(orig_data)
table(orig_data[,'group1'])

ggplot(orig_data[s1], aes(Window, value, group = Promoter)) +
  facet_wrap(~group1, ncol = 2, scales = "free_y") +
  geom_line(color = "grey10", alpha = 0.25, size = 0.1) +
  labs(x = "windows [10bp]", y = "coverage (normalised)") +
  theme_bw()

mean_data <- aggregate(orig_data[, 'value'], list(orig_data$Window, orig_data$group1), mean)
ggplot(mean_data, aes(Group.1, value)) +
  facet_wrap(~Group.2, ncol = 1, scales = "free_y") +
  geom_line(color = "grey10", alpha = 1, size = 1) +
  labs(x = "windows [10bp]", y = "coverage (normalised)") +
  theme_bw()

# clustering
numCores <- parallel::detectCores()

k <- 2:10
t1 <- TSrepr::repr_matrix(matrix, func = repr_sma, args = list(order = 100), normalise = TRUE, func_norm = norm_z)
t1.feaclip <- TSrepr::repr_matrix(matrix, func = repr_feaclip, normalise = TRUE, func_norm = norm_z)
t1.dct <- TSrepr::repr_matrix(matrix, func = repr_dct, args = list(coef = 100), normalise = TRUE, func_norm = norm_z)

clusterings <- parallel::mclapply(k, function(x){p = cluster::pam(t1, x); return(p)}, mc.cores = numCores, mc.cleanup=T)
save(clusterings, file = 'clusterings.rda')

clusterings.feaclip <- parallel::mclapply(k, function(x){p = cluster::pam(t1.feaclip, x); return(p)}, mc.cores = numCores, mc.cleanup=T)
save(clusterings.feaclip, file = 'clusterings.feaclip.rda')

clusterings.dct <- parallel::mclapply(k, function(x){p = cluster::pam(t1.dct, x); return(p)}, mc.cores = numCores, mc.cleanup=T)
save(clusterings.dct, file = 'clusterings.dct.rda')





# inspect clustering results to choose optimal k/N
DB_values <- parallel::mclapply(seq_along(1:length(clusterings)), function(x) {
    ic = clusterCrit::intCriteria(t1, as.integer(clusterings[[x]]$clustering), c("Davies_Bouldin"))
    return(ic)}, 
    mc.cores = numCores, mc.cleanup = TRUE)
save(DB_values, file = 'DB_values.rda')

DB_values.feaclip <- parallel::mclapply(seq_along(1:length(clusterings.feaclip)), function(x) {
  ic = clusterCrit::intCriteria(t1.feaclip, as.integer(clusterings.feaclip[[x]]$clustering), c("Davies_Bouldin"))
  return(ic)}, 
  mc.cores = numCores, mc.cleanup = TRUE)
save(DB_values.feaclip, file = 'DB_values.feaclip.rda')

DB_values.dct <- parallel::mclapply(seq_along(1:length(clusterings.dct)), function(x) {
  ic = clusterCrit::intCriteria(t1.dct, as.integer(clusterings.dct[[x]]$clustering), c("Davies_Bouldin"))
  return(ic)}, 
  mc.cores = numCores, mc.cleanup = TRUE)
save(DB_values.dct, file = 'DB_values.dct.rda')


g1 <- ggplot2::ggplot(data.table(Clusters = k, DBindex = unlist(DB_values)),
       aes(Clusters, DBindex)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  theme_bw()
ggplot2::ggsave(g1, filename = 'g1.pdf')

g1.feaclip <- ggplot2::ggplot(data.table(Clusters = k, DBindex = unlist(DB_values.feaclip)),
                      aes(Clusters, DBindex)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  theme_bw()
ggplot2::ggsave(g1.feaclip, filename = 'g1.feaclip.pdf')

g1.dct <- ggplot2::ggplot(data.table(Clusters = k, DBindex = unlist(DB_values.dct)),
                              aes(Clusters, DBindex)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  theme_bw()
ggplot2::ggsave(g1.dct, filename = 'g1.dct.pdf')

setkey(data_plot, Promoter)
g2 <- ggplot(data_plot[s1], aes(Window, value, group = Promoter)) +
  facet_wrap(~class, ncol = 2, scales = "free_y") +
  geom_line(color = "grey10", alpha = 0.25, size = 0.1) +
  geom_line(data = centers, aes(Window, value),
            color = "firebrick1", alpha = 0.80, size = 1.2) +
  labs(x = "windows [20bp]", y = "coverage (normalised)") +
  theme_bw()
ggsave(g2, filename='g2.pdf')

#feaclip
data_plot <- data.table(melt(data.table(class = as.factor(clusterings[[3]]$clustering),t1.feaclip)))

data_plot[, Window := rep(1:ncol(t1), each = nrow(t1))]
data_plot[, Promoter := rep(mcf10_wt_input_total$computeMatrixRows, ncol(t1))]
centers <- data.table(melt(clusterings[[3]]$medoids))
setnames(centers, c("Var1", "Var2"), c("class", "Window"))
centers[, Promoter := class]


g2.feaclip <- ggplot(data_plot[s1], aes(Window, value, group = Promoter)) +
  facet_wrap(~class, ncol = 2, scales = "free_y") +
  geom_line(color = "grey10", alpha = 0.25, size = 0.1) +
  geom_line(data = centers, aes(Window, value),
            color = "firebrick1", alpha = 0.80, size = 1.2) +
  labs(x = "windows [20bp]", y = "coverage (normalised)") +
  theme_bw()



# attempting the use of seasonal profiling for data reduction


clustering.seas_prof <- parallel::mclapply(k, function(x){p = cluster::pam(matrix_seasprof, x); return(p)}, mc.cores = numCores, mc.cleanup=T)
save(clustering.seas_prof, file = file.path(dest_dir, 'clustering.seas_prof.rda'))

DB_values.seas_prof <- parallel::mclapply(seq_along(1:length(clustering.seas_prof)), function(x) {
  ic = clusterCrit::intCriteria(matrix_seasprof, as.integer(clustering.seas_prof[[x]]$clustering), c("Davies_Bouldin"))
  return(ic)}, 
  mc.cores = numCores, mc.cleanup = TRUE)
save(DB_values.seas_prof, file = file.path(dest_dir, 'DB_values.seas_prof.rda'))

g1.seas_prof <- ggplot2::ggplot(data.table(Clusters = k, DBindex = unlist(DB_values.seas_prof)),
                                aes(Clusters, DBindex)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  theme_bw()
ggplot2::ggsave(g1.seas_prof, filename = file.path(dest_dir, 'g1.seas_prof.pdf'))

selected_cluster <- which.min(unlist(DB_values.seas_prof))

# force 7 from k-manes 
selected_cluster <- 6

data_plot <- data.table(melt(data.table(class = as.factor(clustering.seas_prof[[selected_cluster]]$clustering), matrix_seasprof)))
data_plot[, Window := rep(1:ncol(matrix_seasprof), each = nrow(matrix_seasprof))]
data_plot[, Promoter := rep(mcf10_wt_input_total$computeMatrixRows, ncol(matrix_seasprof))]
centers <- data.table(melt(clustering.seas_prof[[selected_cluster]]$medoids))
centers$Var1 <- rep(c(1:7), ncol(matrix_seasprof))
setnames(centers, c("Var1", "Var2"), c("class", "Window"))
centers[, Promoter := class]

g2.seasprof <- ggplot(data_plot, aes(Window, value)) +
  facet_wrap(~class, ncol = 2, scales = "free_y") +
  geom_line(color = "grey10", alpha = 0.25, size = 0.1) +
  geom_line(data = centers, aes(Window, value),
            color = "firebrick1", alpha = 0.80, size = 1.2) +
  geom_line(data = mean_data, aes(Group.1, value),
            color = "green", size = 1.1) +
  labs(x = "windows [50bp]", y = "coverage (normalised)") +
  theme_bw()

ggplot2::ggsave(g2.seasprof, filename = file.path(dest_dir, 'g2.seasprof.pdf'))

g3 <- ggplot(mean_data, aes(Window, value)) + facet_wrap(~class, ncol = 2, scales = 'free_y') + geom_line(color = 'green')

ggplot(mean_data, aes(Window, value)) + facet_wrap(~class, ncol = 2, scales = 'free_y') + geom_line(color = 'green') + 
  geom_line(data = centers, aes(Window, value),
            color = "firebrick1", alpha = 0.80, size = 1.2)
