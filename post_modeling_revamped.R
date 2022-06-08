library(ridge)
library(glmnet)
library(Matrix)
library(ggplot2)
library(ggrepel)
library(ggnewscale)
library(patchwork)
library(RColorBrewer)
library(LimitBreaks)
plain.theme <- theme(axis.line.x = element_line(color = '#000000', size = 0.5, lineend = 'round'),
                     axis.line.y = element_line(color = '#000000', size = 0.5, lineend = 'round'),
                     axis.ticks = element_line(color = 'black', size = 0.5, lineend = 'round'),
                     axis.ticks.length = unit(4,'points'),
                     panel.background = element_rect(fill = 'white'),
                     text = element_text(size = 12, color = 'black'),
                     axis.text = element_text(size = 12, color = 'black'))
fpkm.counts <- readRDS('rna_fpkm.rds')
normcounts <- readRDS('rna_normcounts.rds')
rna.results <- readRDS('rna_results.rds')
atac.results <- readRDS('atac_results.rds')
pbm <- readRDS('./modeling/pbm.rds')
numtfs <- ncol(pbm)
dir.create(paste0('Pres_plots/'), showWarnings = FALSE)
toexclude <- c('Ighv','Ighj','Ighd','Igkv','Igkj','Iglv','Iglj','Trav','Traj','Trbv','Trbj','Trbd','Trgv','Trgj','Trdv','Trdj','Trdd')
toexclude <- unique(c(unlist(lapply(toexclude,function(i){
  rownames(rna.results)[grep(i,rownames(rna.results))]
}))))
toexclude <- unique(c(toexclude,rownames(fpkm.counts)[which(rowMeans(fpkm.counts) <= mean(fpkm.counts['Cd8a',]))]))
toexclude <- unique(c(toexclude,rownames(fpkm.counts)[which(apply(fpkm.counts,1,median) == 0)]))
rna.results <- rna.results[-match(toexclude,rownames(rna.results)),]
normcounts <- normcounts[-match(toexclude,rownames(normcounts)),]
fpkm.counts <- fpkm.counts[-match(toexclude,rownames(fpkm.counts)),]
sig.genes <- unique(c(which(rna.results$padj_pvn < 0.05),which(rna.results$padj_svn < 0.05),which(rna.results$padj_svp < 0.05)))
z.scorecounts <- do.call(rbind,lapply(seq_along(rownames(normcounts)),function(i){
  meanr <- mean(normcounts[i,])
  sdr <- sd(normcounts[i,])
  if(sdr == 0){
    newrow <- rep(0,ncol(normcounts))
  }
  else{
    newrow <- (normcounts[i,] - meanr)/sdr}
}))
rownames(z.scorecounts) <- rownames(normcounts)
colnames(z.scorecounts) <- colnames(normcounts)
rna.forclust <- z.scorecounts[sig.genes,]
rna.clustcol <- hclust(dist(t(rna.forclust)))
#rna.kmeans <- kmeans(rna.forclust, centers = 7, nstart = 50)
#saveRDS(rna.kmeans, 'rna_kmeans.rds')
rna.kmeans <- readRDS('rna_kmeans.rds')
rna.clustersizes <- table(rna.kmeans$cluster)
rna.sizes.plot <- rna.clustersizes[1]
for(n in 2:length(rna.clustersizes)){rna.sizes.plot <- c(rna.sizes.plot,sum(rna.clustersizes[1:n]))}
rna.clustersizes <- data.frame(rna.clustersizes)
colnames(rna.clustersizes) <- c('Cluster', 'N_genes')
rna.gene.order <- data.frame(Gene = rownames(rna.forclust)[order(rna.kmeans$cluster)], Cluster = rna.kmeans$cluster[order(rna.kmeans$cluster)])
rna.samp.order <- labels(as.dendrogram(rna.clustcol))
rna.kmeans.data <- data.frame(Gene = rownames(rna.forclust), Sample = rep(colnames(rna.forclust), each = nrow(rna.forclust)), Value = as.vector(rna.forclust))
rna.kmeans.plot <- ggplot(data = rna.kmeans.data, aes(Sample, Gene)) +
  geom_tile(aes(fill = Value)) +
  geom_hline(size = 1, yintercept = rna.sizes.plot) +
  scale_fill_gradientn(colours = rev(brewer.pal(9,'RdYlBu'))) +
  scale_x_discrete(limits = c('RNA_G1_neg','RNA_G2_neg','RNA_G3_neg','RNA_G1_pos','RNA_G2_pos','RNA_G3_pos','RNA_G1_stable','RNA_G2_stable','RNA_G3_stable')) +
  scale_y_discrete(limits = rna.gene.order$Gene, labels = NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  plain.theme
#ggsave('./Pres_plots/kmeans_heatmap.pdf', width = 5, height = 8, dpi = 300)
atac.results$RNA.cluster <- unlist(lapply(seq_along(rownames(atac.results)),function(i){
  if(atac.results$Gene[i] %in% rna.gene.order$Gene){paste0('Cluster_', rna.gene.order$Cluster[which(rna.gene.order == atac.results$Gene[i])])}
  else{paste0('No_cluster')}
}))
rna.results$RNA.cluster <- unlist(lapply(seq_along(rownames(rna.results)),function(i){
  if(rna.results$Gene[i] %in% rna.gene.order$Gene){paste0('Cluster_', rna.gene.order$Cluster[which(rna.gene.order == rna.results$Gene[i])])}
  else{paste0('No_cluster')}
}))
atac.norm <- readRDS('atac_normcounts.rds')
atac.z.score <- do.call(rbind,lapply(seq_along(rownames(atac.norm)),function(i){
  meanr <- mean(atac.norm[i,])
  sdr <- sd(atac.norm[i,])
  if(sdr == 0){
    newrow <- rep(0,ncol(atac.norm))
  }
  else{
    newrow <- (atac.norm[i,] - meanr)/sdr}
}))
rownames(atac.z.score) <- rownames(atac.norm)
colnames(atac.z.score) <- colnames(atac.norm)
tdf <- atac.results[which(atac.results$RNA.cluster %in% c("Cluster_6","Cluster_4","Cluster_3","Cluster_2","Cluster_1","Cluster_5","Cluster_7")),]
atac.forclust <- atac.z.score[rownames(tdf),]
atac.peak.order <- data.frame(Peak = rownames(tdf)[order(tdf$RNA.cluster)], Cluster = tdf$RNA.cluster[order(tdf$RNA.cluster)])
atac.kmeans.data <- data.frame(Peak = rownames(atac.forclust), Sample = rep(colnames(atac.forclust), each = nrow(atac.forclust)), Value = as.vector(atac.forclust))
atac.clustersizes <- table(tdf$RNA.cluster)
atac.sizes.plot <- cumsum(atac.clustersizes)
atac.kmeans.plot <- ggplot(data = atac.kmeans.data, aes(Sample, Peak)) +
  geom_tile(aes(fill = Value)) +
  geom_hline(size = 0.5, yintercept = atac.sizes.plot) +
  scale_fill_gradientn(colours = rev(brewer.pal(9,'RdYlBu'))) +
  scale_x_discrete(limits = c('ATAC_G1_neg','ATAC_G2_neg','ATAC_G3_neg','ATAC_G1_pos','ATAC_G2_pos','ATAC_G3_pos','ATAC_G1_stable','ATAC_G2_stable','ATAC_G3_stable')) +
  scale_y_discrete(limits = atac.peak.order$Peak, labels = NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  plain.theme
atac.fpkm <- readRDS('atac_fpkm.rds')
atac.fpkm <- log2(atac.fpkm + 1)
atac.fpkm <- atac.fpkm[rownames(tdf),]
atac.df <- do.call(cbind,lapply(1:ncol(atac.fpkm),function(i){
  tapply(atac.fpkm[,i],tdf$RNA.cluster, mean, na.rm = TRUE)
}))
colnames(atac.df) <- colnames(atac.fpkm)
atac.df <- data.frame(Sample = colnames(atac.df),
                      scale(t(atac.df)))
atac.df <- rbind(atac.df,c('Cor_svp',unlist(sapply(sort(unique(tdf$RNA.cluster)),function(i){
  atac.results.v2 <- data.frame(do.call(cbind,list(tapply(atac.results$lg2FC_svn, atac.results$Gene, mean),tapply(atac.results$lg2FC_svp, atac.results$Gene, mean))))
  colnames(atac.results.v2) <- c('lg2FC_svn','lg2FC_svp')
  atac.results.v2$Gene <- rownames(atac.results.v2)
  both.res <- merge(rna.results[,c(1,5,7,11)], atac.results.v2, by = 'Gene', sort = FALSE, suffixes = c('.rna','.atac'))
  cor(subset(both.res, RNA.cluster == as.character(i))$lg2FC_svp.atac,subset(both.res, RNA.cluster == as.character(i))$lg2FC_svp.rna)
  #summary(lm(lg2FC_svp.atac~lg2FC_svp.rna, data = subset(both.res, RNA.cluster == as.character(i))))$adj.r.squared
}))))
atac.df <- rbind(atac.df,c('Cor_svn',unlist(sapply(sort(unique(tdf$RNA.cluster)),function(i){
  atac.results.v2 <- data.frame(do.call(cbind,list(tapply(atac.results$lg2FC_svn, atac.results$Gene, mean),tapply(atac.results$lg2FC_svp, atac.results$Gene, mean))))
  colnames(atac.results.v2) <- c('lg2FC_svn','lg2FC_svp')
  atac.results.v2$Gene <- rownames(atac.results.v2)
  both.res <- merge(rna.results[,c(1,5,7,11)], atac.results.v2, by = 'Gene', sort = FALSE, suffixes = c('.rna','.atac'))
  cor(subset(both.res, RNA.cluster == as.character(i))$lg2FC_svn.atac,subset(both.res, RNA.cluster == as.character(i))$lg2FC_svn.rna)
  #summary(lm(lg2FC_svn.atac~lg2FC_svn.rna, data = subset(both.res, RNA.cluster == as.character(i))))$adj.r.squared
}))))
svp.cor <- unlist(sapply(1,function(i){
  atac.results.v2 <- data.frame(do.call(cbind,list(tapply(atac.results$lg2FC_svn, atac.results$Gene, mean),tapply(atac.results$lg2FC_svp, atac.results$Gene, mean))))
  colnames(atac.results.v2) <- c('lg2FC_svn','lg2FC_svp')
  atac.results.v2$Gene <- rownames(atac.results.v2)
  both.res <- merge(rna.results[,c(1,5,7,11)], atac.results.v2, by = 'Gene', sort = FALSE, suffixes = c('.rna','.atac'))
  cor(both.res$lg2FC_svp.atac,both.res$lg2FC_svp.rna)
  #summary(lm(lg2FC_svp.atac~lg2FC_svp.rna, data = both.res))$adj.r.squared
}))
svn.cor <- unlist(sapply(1,function(i){
  atac.results.v2 <- data.frame(do.call(cbind,list(tapply(atac.results$lg2FC_svn, atac.results$Gene, mean),tapply(atac.results$lg2FC_svp, atac.results$Gene, mean))))
  colnames(atac.results.v2) <- c('lg2FC_svn','lg2FC_svp')
  atac.results.v2$Gene <- rownames(atac.results.v2)
  both.res <- merge(rna.results[,c(1,5,7,11)], atac.results.v2, by = 'Gene', sort = FALSE, suffixes = c('.rna','.atac'))
  cor(both.res$lg2FC_svn.atac,both.res$lg2FC_svn.rna)
  #summary(lm(lg2FC_svn.atac~lg2FC_svn.rna, data = both.res))$adj.r.squared
}))
atac.df <- data.frame(Celltype = rep(atac.df$Sample,7), bulk_cluster = rep(colnames(atac.df)[2:8], each = nrow(atac.df)), expr = unlist(atac.df[,2:8]))
atac.df$expr <- as.numeric(atac.df$expr)
svp.cor.min <- svp.cor - 1.001*max(abs(min(subset(atac.df, Celltype %in% c('Cor_svp'))$expr)-svp.cor),abs(max(subset(atac.df, Celltype %in% c('Cor_svp'))$expr)-svp.cor))
svp.cor.max <- svp.cor + 1.001*max(abs(min(subset(atac.df, Celltype %in% c('Cor_svp'))$expr)-svp.cor),abs(max(subset(atac.df, Celltype %in% c('Cor_svp'))$expr)-svp.cor))
svn.cor.min <- svn.cor - 1.001*max(abs(min(subset(atac.df, Celltype %in% c('Cor_svn'))$expr)-svn.cor),abs(max(subset(atac.df, Celltype %in% c('Cor_svn'))$expr)-svn.cor))
svn.cor.max <- svn.cor + 1.001*max(abs(min(subset(atac.df, Celltype %in% c('Cor_svn'))$expr)-svn.cor),abs(max(subset(atac.df, Celltype %in% c('Cor_svn'))$expr)-svn.cor))
atac.b.plot <- ggplot(data = subset(atac.df, Celltype %in% colnames(atac.fpkm)), aes(Celltype, bulk_cluster)) +
  geom_tile(aes(fill = expr), color = "#FFFFFF", size = 0.3) +
  scale_fill_gradientn(colors = rev(brewer.pal(11,'RdYlBu')), limits = c(-3,3)) +
  new_scale_fill() +
  geom_tile(data = subset(atac.df, Celltype %in% c('Cor_svp')), aes(x = Celltype, y = bulk_cluster, fill = expr), color = "#FFFFFF", size = 0.3) +
  scale_fill_gradient2(high = brewer.pal(11,'RdBu')[1], mid = brewer.pal(11,'RdBu')[6], low = brewer.pal(11,'RdBu')[11], limits = c(svp.cor.min,svp.cor.max), midpoint = svp.cor) +
  new_scale_fill() +
  geom_tile(data = subset(atac.df, Celltype %in% c('Cor_svn')), aes(x = Celltype, y = bulk_cluster, fill = expr), color = "#FFFFFF", size = 0.3) +
  scale_fill_gradient2(high = brewer.pal(11,'RdBu')[1], mid = brewer.pal(11,'RdBu')[6], low = brewer.pal(11,'RdBu')[11], limits = c(svn.cor.min,svn.cor.max), midpoint = svn.cor) +
  scale_x_discrete(limits = c('ATAC_G1_neg','ATAC_G2_neg','ATAC_G3_neg','ATAC_G1_pos','ATAC_G2_pos','ATAC_G3_pos','ATAC_G1_stable','ATAC_G2_stable','ATAC_G3_stable',NA,'Cor_svp','Cor_svn')) +
  scale_y_discrete(limits = sort(unique(tdf$RNA.cluster))) + 
  plain.theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
rna.fpkm <- log2(fpkm.counts + 1)
rna.fpkm <- rna.fpkm[names(rna.kmeans$cluster),]
rna.df <- do.call(cbind,lapply(1:ncol(rna.fpkm),function(i){
  tapply(rna.fpkm[,i],rna.kmeans$cluster, mean, na.rm = TRUE)
}))
rownames(rna.df) <- paste("Cluster",as.vector(tapply(rna.kmeans$cluster,rna.kmeans$cluster, unique)), sep = "_")
colnames(rna.df) <- colnames(rna.fpkm)
rna.df <- data.frame(Sample = colnames(rna.df),
                      scale(t(rna.df)))
rna.df <- data.frame(Celltype = rep(rna.df$Sample,7), bulk_cluster = rep(colnames(rna.df)[2:8], each = nrow(rna.df)), expr = unlist(rna.df[,2:8]))
rna.b.plot <- ggplot(data = rna.df, aes(Celltype, bulk_cluster)) +
  geom_tile(aes(fill = expr), color = "#FFFFFF", size = 0.3) +
  scale_fill_gradientn(colors = rev(brewer.pal(11,'RdYlBu')), limits = c(-3,3)) +
  scale_x_discrete(limits = c('RNA_G1_neg','RNA_G2_neg','RNA_G3_neg','RNA_G1_pos','RNA_G2_pos','RNA_G3_pos','RNA_G1_stable','RNA_G2_stable','RNA_G3_stable')) +
  plain.theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
pra <- rna.b.plot + atac.b.plot + patchwork::plot_layout(ncol = 2, widths = c(1,1.35))
ggsave('./Pres_plots/atac_as_bulk_new.pdf', height = 5, width = 12, dpi = 300)

rna.clustersizes$Cluster <- paste0('Cluster_',1:7)
rna.clustersizes$neg_expr <- unlist(lapply(1:nrow(rna.clustersizes),function(i){
  mean(rna.forclust[which(rna.kmeans$cluster == i),grep('neg', colnames(rna.forclust))])
}))
rna.clustersizes$pos_expr <- unlist(lapply(1:nrow(rna.clustersizes),function(i){
  mean(rna.forclust[which(rna.kmeans$cluster == i),grep('pos', colnames(rna.forclust))])
}))
rna.clustersizes$stable_expr <- unlist(lapply(1:nrow(rna.clustersizes),function(i){
  mean(rna.forclust[which(rna.kmeans$cluster == i),grep('stable', colnames(rna.forclust))])
}))
rna.clustersizes$designation <- rep(NA,nrow(rna.clustersizes))
rna.clustersizes$designation[which(rna.clustersizes$neg_expr > rna.clustersizes$pos_expr & rna.clustersizes$pos_expr > rna.clustersizes$stable_expr)] <- 'decreasing'
rna.clustersizes$designation[which.max(rna.clustersizes$stable_expr)] <- 'stable_high'
rna.clustersizes$designation[setdiff(which(rna.clustersizes$stable_expr > rna.clustersizes$pos_expr & rna.clustersizes$pos_expr > rna.clustersizes$neg_expr),which.max(rna.clustersizes$stable_expr))] <- 'increasing'
rna.clustersizes$designation[which.max(rna.clustersizes$pos_expr)] <- 'pos_high'
rna.clustersizes$designation[which.max(rna.clustersizes$neg_expr)] <- 'neg_high'
rna.clustersizes$designation[setdiff(which(rna.clustersizes$stable_expr > rna.clustersizes$pos_expr & rna.clustersizes$neg_expr > rna.clustersizes$pos_expr),which.max(rna.clustersizes$neg_expr))] <- 'pos_low'
rna.clustersizes$designation[setdiff(which(rna.clustersizes$pos_expr > rna.clustersizes$stable_expr & rna.clustersizes$pos_expr > rna.clustersizes$neg_expr),which.max(rna.clustersizes$pos_expr))] <- 'pos_peak'
rna.clustersizes <- rbind(data.frame(Cluster = 'all',
                                     N_genes = length(unique(atac.results$Gene)),
                                     neg_expr = mean(z.scorecounts[which(rownames(z.scorecounts) %in% unique(atac.results$Gene)),grep('neg', colnames(z.scorecounts))]),
                                     pos_expr = mean(z.scorecounts[which(rownames(z.scorecounts) %in% unique(atac.results$Gene)),grep('pos', colnames(z.scorecounts))]),
                                     stable_expr = mean(z.scorecounts[which(rownames(z.scorecounts) %in% unique(atac.results$Gene)),grep('stable', colnames(z.scorecounts))]),
                                     designation = 'all_peaks',
                                     stringsAsFactors = FALSE),
                          rna.clustersizes)
clust.dec <- rna.clustersizes$Cluster[which(rna.clustersizes$designation == 'decreasing')]
clust.inc <- rna.clustersizes$Cluster[which(rna.clustersizes$designation == 'increasing')]
#saveRDS(rna.clustersizes,'rna_cluster_info.rds')
famkey <- readRDS(paste0('./modeling/famkey.rds'))
famkey$fam_members <- as.character(famkey$fam_members)
svn.model <- readRDS('./modeling/svn_ridge_model.rds')
svp.model <- readRDS('./modeling/svp_ridge_model.rds')
svn.pvals <- readRDS('./modeling/svn_ridge_pvals.rds')
svp.pvals <- readRDS('./modeling/svp_ridge_pvals.rds')
TF.coeffs <- data.frame(TFs = names(coef(svn.model)[2:(1+numtfs)]), stringsAsFactors = FALSE)
TF.coeffs$svn_coefs <- coef(svn.model)[2:(1+numtfs)]
TF.coeffs$svp_coefs <- coef(svp.model)[2:(1+numtfs)]
TF.coeffs$svn_pvals <- svn.pvals$pval[,svn.pvals$chosen.nPCs]
TF.coeffs$svp_pvals <- svp.pvals$pval[,svp.pvals$chosen.nPCs]
sig.res <- subset(TF.coeffs, svn_pvals < 0.05 | svp_pvals < 0.05)
rownames(sig.res) <- sig.res$TFs
rna.clustersizes$svp_mse <- unlist(lapply(1:nrow(rna.clustersizes), function(h){
  currclust <- as.character(rna.clustersizes$Cluster[h])
  if(currclust == 'all'){clusterpeaks <- atac.results$peak}
  else{clusterpeaks <- atac.results$peak[which(atac.results$RNA.cluster == currclust)]}
  clusterpeaks <- intersect(clusterpeaks, rownames(pbm))
  svp_var_test <- atac.results[clusterpeaks,'lg2FC_svp']
  x_vars_test <- pbm[clusterpeaks,]
  x_vars_test <- predict(svp.model, newdata = data.frame(x_vars_test))
  mean((svp_var_test - x_vars_test)^2)
}))
rna.clustersizes$svn_mse <- unlist(lapply(1:nrow(rna.clustersizes), function(h){
  currclust <- as.character(rna.clustersizes$Cluster[h])
  if(currclust == 'all'){clusterpeaks <- atac.results$peak}
  else{clusterpeaks <- atac.results$peak[which(atac.results$RNA.cluster == currclust)]}
  clusterpeaks <- intersect(clusterpeaks, rownames(pbm))
  svn_var_test <- atac.results[clusterpeaks,'lg2FC_svn']
  x_vars_test <- pbm[clusterpeaks,]
  x_vars_test <- predict(svn.model, newdata = data.frame(x_vars_test))
  mean((svn_var_test - x_vars_test)^2)
}))
svp.models <- readRDS(paste0('./modeling/svp_ridge_test_models.rds'))
svn.models <- readRDS(paste0('./modeling/svn_ridge_test_models.rds'))
pbm.tf <- do.call(cbind,lapply(1:ncol(pbm),function(i){
  pbm[,i] > 0
}))
colnames(pbm.tf) <- colnames(pbm)
famkey <- cbind(famkey, do.call(rbind,lapply(1:nrow(famkey),function(i){
  cf <- famkey$fam_id[i]
  table(atac.results$peak_type[which(atac.results$peak %in% names(which(pbm.tf[,cf])))])/length(which(atac.results$peak %in% names(which(pbm.tf[,cf]))))
})))
tt_fam <- strsplit(famkey$fam_members, split = ",")
names(tt_fam) <- famkey$fam_id
tfs.in.model <- unique(unlist(strsplit(famkey$fam_members, split = ",")))
system.time(sig.res <- cbind(sig.res, do.call(rbind,lapply(1:nrow(sig.res),function(i){
  currtf <-  sig.res$TFs[i]
  currtfn <- which(colnames(pbm) == currtf)
  modelmod <- svp.models[[currtf]]
  infs <- do.call(cbind,lapply(1:nrow(rna.clustersizes), function(h){
    currclust <- as.character(rna.clustersizes$Cluster[h])
    if(currclust == 'all'){clusterpeaks <- atac.results$peak}
    else{clusterpeaks <- atac.results$peak[which(atac.results$RNA.cluster == currclust)]}
    clusterpeaks <- intersect(clusterpeaks, rownames(pbm))
    if(rna.clustersizes[h,5] - rna.clustersizes[h,4] > 0){clusterpeaks <- intersect(clusterpeaks,atac.results$peak[which(atac.results$lg2FC_svp > 0)])}
    else{clusterpeaks <- intersect(clusterpeaks,atac.results$peak[which(atac.results$lg2FC_svp < 0)])}
    svp_var_test <- atac.results[clusterpeaks,'lg2FC_svp']
    x_vars_test <- pbm[clusterpeaks,]
    x_vars_test <- predict(svp.model, newdata = data.frame(x_vars_test))
    cor.o <- cor(svp_var_test, x_vars_test)
    x_vars_test <- pbm[clusterpeaks,-currtfn]
    x_vars_test <- predict(modelmod, newdata = data.frame(x_vars_test))
    cor.n <- cor(svp_var_test, x_vars_test)
    cor.change <- cor.o - cor.n
  }))
  colnames(infs) <- paste0(rna.clustersizes$Cluster,'_svp')
  infs
}))))
system.time(sig.res <- cbind(sig.res, do.call(rbind,lapply(1:nrow(sig.res),function(i){
  currtf <-  sig.res$TFs[i]
  currtfn <- which(colnames(pbm) == currtf)
  modelmod <- svn.models[[currtf]]
  infs <- do.call(cbind,lapply(1:nrow(rna.clustersizes), function(h){
    currclust <- as.character(rna.clustersizes$Cluster[h])
    if(currclust == 'all'){clusterpeaks <- atac.results$peak}
    else{clusterpeaks <- atac.results$peak[which(atac.results$RNA.cluster == currclust)]}
    clusterpeaks <- intersect(clusterpeaks, rownames(pbm))
    if(rna.clustersizes[h,5] - rna.clustersizes[h,3] > 0){clusterpeaks <- intersect(clusterpeaks,atac.results$peak[which(atac.results$lg2FC_svn > 0)])}
    else{clusterpeaks <- intersect(clusterpeaks,atac.results$peak[which(atac.results$lg2FC_svn < 0)])}
    svn_var_test <- atac.results[clusterpeaks,'lg2FC_svn']
    x_vars_test <- pbm[clusterpeaks,]
    x_vars_test <- predict(svn.model, newdata = data.frame(x_vars_test))
    cor.o <- cor(svn_var_test, x_vars_test)
    x_vars_test <- pbm[clusterpeaks,-currtfn]
    x_vars_test <- predict(modelmod, newdata = data.frame(x_vars_test))
    cor.n <- cor(svn_var_test, x_vars_test)
    cor.change <- cor.o - cor.n
  }))
  colnames(infs) <- paste0(rna.clustersizes$Cluster,'_svn')
  infs
}))))
totest <- c('pos_up','neg_up')
system.time(sig.res <- cbind(sig.res, do.call(rbind,lapply(1:nrow(sig.res),function(i){
  currtf <-  sig.res$TFs[i]
  currtfn <- which(colnames(pbm) == currtf)
  modelmod <- svp.models[[currtf]]
  infs <- do.call(cbind,lapply(seq_along(totest), function(h){
    temp1 <- read.delim(paste0('../../Il10_scRNAseq/Final_version/',totest[h],'.txt'), header = FALSE)[,1]
    clusterpeaks <- atac.results$peak[which(atac.results$Gene %in% temp1)]
    clusterpeaks <- intersect(clusterpeaks, rownames(pbm))
    if(mean(rna.results$lg2FC_svp[which(rna.results$Gene %in% temp1)]) > 0){clusterpeaks <- intersect(clusterpeaks,atac.results$peak[which(atac.results$lg2FC_svp > 0)])}
    else{clusterpeaks <- intersect(clusterpeaks,atac.results$peak[which(atac.results$lg2FC_svp < 0)])}
    svp_var_test <- atac.results[clusterpeaks,'lg2FC_svp']
    x_vars_test <- pbm[clusterpeaks,]
    x_vars_test <- predict(svp.model, newdata = data.frame(x_vars_test))
    cor.o <- cor(svp_var_test, x_vars_test)
    x_vars_test <- pbm[clusterpeaks,-currtfn]
    x_vars_test <- predict(modelmod, newdata = data.frame(x_vars_test))
    cor.n <- cor(svp_var_test, x_vars_test)
    cor.change <- cor.o - cor.n
  }))
  colnames(infs) <- paste0(totest,'_svp')
  infs
}))))
system.time(sig.res <- cbind(sig.res, do.call(rbind,lapply(1:nrow(sig.res),function(i){
  currtf <-  sig.res$TFs[i]
  currtfn <- which(colnames(pbm) == currtf)
  modelmod <- svn.models[[currtf]]
  infs <- do.call(cbind,lapply(seq_along(totest), function(h){
    temp1 <- read.delim(paste0('../../Il10_scRNAseq/Final_version/',totest[h],'.txt'), header = FALSE)[,1]
    clusterpeaks <- atac.results$peak[which(atac.results$Gene %in% temp1)]
    clusterpeaks <- intersect(clusterpeaks, rownames(pbm))
    if(mean(rna.results$lg2FC_svn[which(rna.results$Gene %in% temp1)]) > 0){clusterpeaks <- intersect(clusterpeaks,atac.results$peak[which(atac.results$lg2FC_svn > 0)])}
    else{clusterpeaks <- intersect(clusterpeaks,atac.results$peak[which(atac.results$lg2FC_svn < 0)])}
    svn_var_test <- atac.results[clusterpeaks,'lg2FC_svn']
    x_vars_test <- pbm[clusterpeaks,]
    x_vars_test <- predict(svn.model, newdata = data.frame(x_vars_test))
    cor.o <- cor(svn_var_test, x_vars_test)
    x_vars_test <- pbm[clusterpeaks,-currtfn]
    x_vars_test <- predict(modelmod, newdata = data.frame(x_vars_test))
    cor.n <- cor(svn_var_test, x_vars_test)
    cor.change <- cor.o - cor.n
  }))
  colnames(infs) <- paste0(totest,'_svn')
  infs
}))))
secreted <- unique(read.delim('secreted.txt', stringsAsFactors = FALSE, header = TRUE)[,3])
secreted <- secreted[which(nchar(secreted) > 0)]
surface <- unique(read.delim('surface.txt', stringsAsFactors = FALSE, header = TRUE)[,3])
surface <- surface[which(nchar(surface) > 0)]
surface <- setdiff(surface,secreted)
stable.sec <- intersect(rna.results$Gene[which(rna.results$RNA.cluster == 'Cluster_2' & rna.results$padj_svn < 0.05 & rna.results$lg2FC_svn > 0.5)],secreted)
stable.sec <- unique(c(stable.sec,intersect(rna.results$Gene[which(rna.results$RNA.cluster == 'Cluster_4' & rna.results$padj_svn < 0.05 & rna.results$lg2FC_svn > 0.5)],secreted)))
stable.sur <- intersect(rna.results$Gene[which(rna.results$RNA.cluster == 'Cluster_2' & rna.results$padj_svn < 0.05 & rna.results$lg2FC_svn > 0.5)],surface)
stable.sur <- unique(c(stable.sur,intersect(rna.results$Gene[which(rna.results$RNA.cluster == 'Cluster_4' & rna.results$padj_svn < 0.05 & rna.results$lg2FC_svn > 0.5)],surface)))
stable_int <- sort(c(stable.sur,stable.sec))
neg.sec <- intersect(rna.results$Gene[which(rna.results$RNA.cluster == 'Cluster_1' & rna.results$padj_svn < 0.05 & rna.results$lg2FC_svn < -0.5)],secreted)
neg.sec <- unique(c(neg.sec,intersect(rna.results$Gene[which(rna.results$RNA.cluster == 'Cluster_5' & rna.results$padj_svn < 0.05 & rna.results$lg2FC_svn < -0.5)],secreted)))
neg.sur <- intersect(rna.results$Gene[which(rna.results$RNA.cluster == 'Cluster_1' & rna.results$padj_svn < 0.05 & rna.results$lg2FC_svn < -0.5)],surface)
neg.sur <- unique(c(neg.sur,intersect(rna.results$Gene[which(rna.results$RNA.cluster == 'Cluster_5' & rna.results$padj_svn < 0.05 & rna.results$lg2FC_svn < -0.5)],surface)))
neg_int <- sort(c(neg.sur,neg.sec))
totest2 <- list(stable_int,neg_int)
names(totest2) <- c('stable_int','neg_int')
system.time(sig.res <- cbind(sig.res, do.call(rbind,lapply(1:nrow(sig.res),function(i){
  currtf <-  sig.res$TFs[i]
  currtfn <- which(colnames(pbm) == currtf)
  modelmod <- svp.models[[currtf]]
  infs <- do.call(cbind,lapply(seq_along(totest2), function(h){
    temp1 <- totest2[[h]]
    clusterpeaks <- atac.results$peak[which(atac.results$Gene %in% temp1)]
    clusterpeaks <- intersect(clusterpeaks, rownames(pbm))
    if(mean(rna.results$lg2FC_svp[which(rna.results$Gene %in% temp1)]) > 0){clusterpeaks <- intersect(clusterpeaks,atac.results$peak[which(atac.results$lg2FC_svp > 0)])}
    else{clusterpeaks <- intersect(clusterpeaks,atac.results$peak[which(atac.results$lg2FC_svp < 0)])}
    svp_var_test <- atac.results[clusterpeaks,'lg2FC_svp']
    x_vars_test <- pbm[clusterpeaks,]
    x_vars_test <- predict(svp.model, newdata = data.frame(x_vars_test))
    cor.o <- cor(svp_var_test, x_vars_test)
    x_vars_test <- pbm[clusterpeaks,-currtfn]
    x_vars_test <- predict(modelmod, newdata = data.frame(x_vars_test))
    cor.n <- cor(svp_var_test, x_vars_test)
    cor.change <- cor.o - cor.n
  }))
  colnames(infs) <- paste0(names(totest2),'_svp')
  infs
}))))
system.time(sig.res <- cbind(sig.res, do.call(rbind,lapply(1:nrow(sig.res),function(i){
  currtf <-  sig.res$TFs[i]
  currtfn <- which(colnames(pbm) == currtf)
  modelmod <- svn.models[[currtf]]
  infs <- do.call(cbind,lapply(seq_along(totest2), function(h){
    temp1 <- totest2[[h]]
    clusterpeaks <- atac.results$peak[which(atac.results$Gene %in% temp1)]
    clusterpeaks <- intersect(clusterpeaks, rownames(pbm))
    if(mean(rna.results$lg2FC_svn[which(rna.results$Gene %in% temp1)]) > 0){clusterpeaks <- intersect(clusterpeaks,atac.results$peak[which(atac.results$lg2FC_svn > 0)])}
    else{clusterpeaks <- intersect(clusterpeaks,atac.results$peak[which(atac.results$lg2FC_svn < 0)])}
    svn_var_test <- atac.results[clusterpeaks,'lg2FC_svn']
    x_vars_test <- pbm[clusterpeaks,]
    x_vars_test <- predict(svn.model, newdata = data.frame(x_vars_test))
    cor.o <- cor(svn_var_test, x_vars_test)
    x_vars_test <- pbm[clusterpeaks,-currtfn]
    x_vars_test <- predict(modelmod, newdata = data.frame(x_vars_test))
    cor.n <- cor(svn_var_test, x_vars_test)
    cor.change <- cor.o - cor.n
  }))
  colnames(infs) <- paste0(names(totest2),'_svn')
  infs
}))))
system.time(sig.res <- cbind(sig.res, do.call(rbind,lapply(1:nrow(sig.res),function(i){
  currtf <-  sig.res$TFs[i]
  currtfn <- which(colnames(pbm) == currtf)
  infs <- do.call(cbind,lapply(1:nrow(rna.clustersizes), function(h){
    currclust <- as.character(rna.clustersizes$Cluster[h])
    if(currclust == 'all'){clusterpeaks <- atac.results$peak}
    else{clusterpeaks <- atac.results$peak[which(atac.results$RNA.cluster == currclust)]}
    clusterpeaks <- intersect(clusterpeaks, names(which(pbm[,currtfn] > 0)))
    mean(atac.results[clusterpeaks,]$lg2FC_svp)
  }))
  colnames(infs) <- paste0(rna.clustersizes$Cluster,'_svp_fc')
  infs
}))))
system.time(sig.res <- cbind(sig.res, do.call(rbind,lapply(1:nrow(sig.res),function(i){
  currtf <-  sig.res$TFs[i]
  currtfn <- which(colnames(pbm) == currtf)
  infs <- do.call(cbind,lapply(1:nrow(rna.clustersizes), function(h){
    currclust <- as.character(rna.clustersizes$Cluster[h])
    if(currclust == 'all'){clusterpeaks <- atac.results$peak}
    else{clusterpeaks <- atac.results$peak[which(atac.results$RNA.cluster == currclust)]}
    clusterpeaks <- intersect(clusterpeaks, names(which(pbm[,currtfn] > 0)))
    mean(atac.results[clusterpeaks,]$lg2FC_svn)
  }))
  colnames(infs) <- paste0(rna.clustersizes$Cluster,'_svn_fc')
  infs
}))))
system.time(sig.res <- cbind(sig.res, do.call(rbind,lapply(1:nrow(sig.res),function(i){
  currtf <-  sig.res$TFs[i]
  currtfn <- which(colnames(pbm) == currtf)
  modelmod <- svp.models[[currtf]]
  infs <- do.call(cbind,lapply(seq_along(totest), function(h){
    temp1 <- read.delim(paste0('../../Il10_scRNAseq/Final_version/',totest[h],'.txt'), header = FALSE)[,1]
    clusterpeaks <- atac.results$peak[which(atac.results$Gene %in% temp1)]
    clusterpeaks <- intersect(clusterpeaks, names(which(pbm[,currtfn] > 0)))
    mean(atac.results[clusterpeaks,]$lg2FC_svp)
  }))
  colnames(infs) <- paste0(totest,'_svp_fc')
  infs
}))))
system.time(sig.res <- cbind(sig.res, do.call(rbind,lapply(1:nrow(sig.res),function(i){
  currtf <-  sig.res$TFs[i]
  currtfn <- which(colnames(pbm) == currtf)
  modelmod <- svn.models[[currtf]]
  infs <- do.call(cbind,lapply(seq_along(totest), function(h){
    temp1 <- read.delim(paste0('../../Il10_scRNAseq/Final_version/',totest[h],'.txt'), header = FALSE)[,1]
    clusterpeaks <- atac.results$peak[which(atac.results$Gene %in% temp1)]
    clusterpeaks <- intersect(clusterpeaks, names(which(pbm[,currtfn] > 0)))
    mean(atac.results[clusterpeaks,]$lg2FC_svn)
  }))
  colnames(infs) <- paste0(totest,'_svn_fc')
  infs
}))))
system.time(sig.res <- cbind(sig.res, do.call(rbind,lapply(1:nrow(sig.res),function(i){
  currtf <-  sig.res$TFs[i]
  currtfn <- which(colnames(pbm) == currtf)
  modelmod <- svp.models[[currtf]]
  infs <- do.call(cbind,lapply(seq_along(totest2), function(h){
    temp1 <- totest2[[h]]
    clusterpeaks <- atac.results$peak[which(atac.results$Gene %in% temp1)]
    clusterpeaks <- intersect(clusterpeaks, names(which(pbm[,currtfn] > 0)))
    mean(atac.results[clusterpeaks,]$lg2FC_svp)
  }))
  colnames(infs) <- paste0(names(totest2),'_svp_fc')
  infs
}))))
system.time(sig.res <- cbind(sig.res, do.call(rbind,lapply(1:nrow(sig.res),function(i){
  currtf <-  sig.res$TFs[i]
  currtfn <- which(colnames(pbm) == currtf)
  modelmod <- svn.models[[currtf]]
  infs <- do.call(cbind,lapply(seq_along(totest2), function(h){
    temp1 <- totest2[[h]]
    clusterpeaks <- atac.results$peak[which(atac.results$Gene %in% temp1)]
    clusterpeaks <- intersect(clusterpeaks, names(which(pbm[,currtfn] > 0)))
    mean(atac.results[clusterpeaks,]$lg2FC_svn)
  }))
  colnames(infs) <- paste0(names(totest2),'_svn_fc')
  infs
}))))
plotdata <- merge(famkey,sig.res, by.x = 'fam_id', by.y = 'TFs')
svn.cutoff <- sd(plotdata$svn_coefs)
svp.cutoff <- sd(plotdata$svp_coefs)
pval.min <- min(min(sig.res$svp_pvals[which(sig.res$svp_pvals > 0)]),min(sig.res$svn_pvals[which(sig.res$svn_pvals > 0)]))*.1
plotdata$svn_pvals <- -log10(plotdata$svn_pvals + pval.min)
plotdata$svp_pvals <- -log10(plotdata$svp_pvals + pval.min)
co1 <- subset(plotdata, abs(svn_coefs) > svn.cutoff | abs(svp_coefs) > svp.cutoff)
toplot <- colnames(plotdata)[grep("_svp", colnames(plotdata))]
toplotfc <- toplot[grep('_fc',toplot)]
toplot <- toplot[-grep('_fc',toplot)]
fams.svp <- do.call(rbind,lapply(toplot, function(i){
  currcol <- which(colnames(plotdata) == as.character(i))
  fccol <- which(colnames(plotdata) == paste0(as.character(i),'_fc'))
  tmpdata1 <- plotdata[,c(1:6,currcol,fccol)]
  colnames(tmpdata1) <- c(colnames(plotdata)[1:6],'curr_col','fc_col')
  coval <- tmpdata1[order(tmpdata1$curr_col, decreasing = TRUE),]
  covalt <- 0 + sd(tmpdata1$curr_col)
  if(length(which(coval$curr_col > covalt)) > 3){
    coval <- coval[which(coval$curr_col > covalt),]
    coval <- coval$fam_id}
  else{coval <- coval$fam_id[1:3]}
  coval <- intersect(tmpdata1$fam_id[which(tmpdata1$curr_col > 0)],coval)
  tmpdata1 <- tmpdata1[which(tmpdata1$svp_pvals > -log10(0.01)),]
  tmpdata <- tmpdata1[which(tmpdata1$fam_id %in% coval),]
  tmpdata <- tmpdata[,c(1,2)]
  tmpdata$Cluster <- substring(as.character(i), 1, (regexpr('svp',as.character(i))-2))
  tmpdata$comp <- 'svp'
  tmpdata
}))
toplot.plots.svp <- lapply(toplot, function(i){
  currcol <- which(colnames(plotdata) == as.character(i))
  fccol <- which(colnames(plotdata) == paste0(as.character(i),'_fc'))
  tmpdata1 <- plotdata[,c(1:6,currcol,fccol)]
  colnames(tmpdata1) <- c(colnames(plotdata)[1:6],'curr_col','fc_col')
  coval <- tmpdata1[order(tmpdata1$curr_col, decreasing = TRUE),]
  covalt <- 0 + sd(tmpdata1$curr_col)
  if(length(which(coval$curr_col > covalt)) > 3){
    coval <- coval[which(coval$curr_col > covalt),]
    coval <- coval$fam_id}
  else{coval <- coval$fam_id[1:3]}
  coval <- intersect(tmpdata1$fam_id[which(tmpdata1$curr_col > 0)],coval)
  tmpdata1 <- tmpdata1[which(tmpdata1$svp_pvals > -log10(0.01)),]
  tmpdata <- tmpdata1[which(tmpdata1$fam_id %in% coval),]
  #tmpdata1 <- tmpdata1[which(tmpdata1$curr_col > 0),]
  if(grepl('all', as.character(i))){curmean <- mean(atac.results$lg2FC_svp)}
  else{if(grepl('Cluster', as.character(i))){curmean <- substring(as.character(i),1,gregexpr('_',as.character(i))[[1]][2]-1)
    curmean <- mean(atac.results$lg2FC_svp[which(atac.results$RNA.cluster == curmean)])}
  else{
    if(grepl('int', as.character(i))){curmean <- substring(as.character(i),1,gregexpr('_',as.character(i))[[1]][2]-1)
    curmean <- totest2[[curmean]]
    curmean <- mean(atac.results$lg2FC_svp[which(atac.results$Gene %in% curmean)])}
    else{
    curmean <- substring(as.character(i),1,gregexpr('_',as.character(i))[[1]][2]-1)
    curmean <- read.delim(paste0('../../Il10_scRNAseq/Final_version/',curmean,'.txt'), header = FALSE)[,1]
    curmean <- mean(atac.results$lg2FC_svp[which(atac.results$Gene %in% curmean)])}}}
  xmax <- max(abs(tmpdata1$curr_col))
  ymax <- max(abs(tmpdata1$fc_col - curmean))
  ymin <- curmean - ymax
  ymax <- curmean + ymax
  ggplot(tmpdata1, aes(curr_col, fc_col)) +
    geom_hline(yintercept = curmean, linetype = 2, size = 0.5) +
    geom_vline(xintercept = 0, linetype = 2, size = 0.5) +
    geom_point(aes(fill = svp_coefs), shape = 21, size = 3) +
    scale_fill_gradient2(midpoint = 0, low = brewer.pal(9,'RdYlBu')[9], mid  = brewer.pal(9,'RdYlBu')[5], high = brewer.pal(9,'RdYlBu')[1], na.value = "#FFFFFF", name = as.character(i)) +
    geom_text_repel(data = tmpdata, aes(label = fam_members), size = 4, fontface = "italic", color = "#005500", box.padding = .4, min.segment.length = unit(0.01,'lines'), segment.size = 0.3) +
    xlim(c(-xmax,xmax)) +
    ylim(c(ymin,ymax)) +
    plain.theme})
names(toplot.plots.svp) <- toplot
toplot <- colnames(plotdata)[grep("_svn", colnames(plotdata))]
toplotfc <- toplot[grep('_fc',toplot)]
toplot <- toplot[-grep('_fc',toplot)]
fams.svn <- do.call(rbind,lapply(toplot, function(i){
  currcol <- which(colnames(plotdata) == as.character(i))
  fccol <- which(colnames(plotdata) == paste0(as.character(i),'_fc'))
  tmpdata1 <- plotdata[,c(1:6,currcol,fccol)]
  colnames(tmpdata1) <- c(colnames(plotdata)[1:6],'curr_col','fc_col')
  coval <- tmpdata1[order(tmpdata1$curr_col, decreasing = TRUE),]
  covalt <- 0 + sd(tmpdata1$curr_col)
  if(length(which(coval$curr_col > covalt)) > 3){
    coval <- coval[which(coval$curr_col > covalt),]
    coval <- coval$fam_id}
  else{coval <- coval$fam_id[1:3]}
  coval <- intersect(tmpdata1$fam_id[which(tmpdata1$curr_col > 0)],coval)
  tmpdata1 <- tmpdata1[which(tmpdata1$svn_pvals > -log10(0.01)),]
  tmpdata <- tmpdata1[which(tmpdata1$fam_id %in% coval),]
  tmpdata <- tmpdata[,c(1,2)]
  tmpdata$Cluster <- substring(as.character(i), 1, (regexpr('svn',as.character(i))-2))
  tmpdata$comp <- 'svn'
  tmpdata
  }))
toplot.plots.svn <- lapply(toplot, function(i){
  currcol <- which(colnames(plotdata) == as.character(i))
  fccol <- which(colnames(plotdata) == paste0(as.character(i),'_fc'))
  tmpdata1 <- plotdata[,c(1:6,currcol,fccol)]
  colnames(tmpdata1) <- c(colnames(plotdata)[1:6],'curr_col','fc_col')
  coval <- tmpdata1[order(tmpdata1$curr_col, decreasing = TRUE),]
  covalt <- 0 + sd(tmpdata1$curr_col)
  if(length(which(coval$curr_col > covalt)) > 3){
    coval <- coval[which(coval$curr_col > covalt),]
    coval <- coval$fam_id}
  else{coval <- coval$fam_id[1:3]}
  coval <- intersect(tmpdata1$fam_id[which(tmpdata1$curr_col > 0)],coval)
  tmpdata1 <- tmpdata1[which(tmpdata1$svn_pvals > -log10(0.01)),]
  tmpdata <- tmpdata1[which(tmpdata1$fam_id %in% coval),]
  #tmpdata1 <- tmpdata1[which(tmpdata1$curr_col > 0),]
  if(grepl('all', as.character(i))){curmean <- mean(atac.results$lg2FC_svn)}
  else{if(grepl('Cluster', as.character(i))){curmean <- substring(as.character(i),1,gregexpr('_',as.character(i))[[1]][2]-1)
  curmean <- mean(atac.results$lg2FC_svn[which(atac.results$RNA.cluster == curmean)])}
    else{
      if(grepl('int', as.character(i))){curmean <- substring(as.character(i),1,gregexpr('_',as.character(i))[[1]][2]-1)
      curmean <- totest2[[curmean]]
      curmean <- mean(atac.results$lg2FC_svn[which(atac.results$Gene %in% curmean)])}
      else{
        curmean <- substring(as.character(i),1,gregexpr('_',as.character(i))[[1]][2]-1)
        curmean <- read.delim(paste0('../../Il10_scRNAseq/Final_version/',curmean,'.txt'), header = FALSE)[,1]
        curmean <- mean(atac.results$lg2FC_svn[which(atac.results$Gene %in% curmean)])}}}
  xmax <- max(abs(tmpdata1$curr_col))
  ymax <- max(abs(tmpdata1$fc_col - curmean))
  ymin <- curmean - ymax
  ymax <- curmean + ymax
  ggplot(tmpdata1, aes(curr_col, fc_col)) +
    geom_hline(yintercept = curmean, linetype = 2, size = 0.5) +
    geom_vline(xintercept = 0, linetype = 2, size = 0.5) +
    geom_point(aes(fill = svn_coefs), shape = 21, size = 3) +
    scale_fill_gradient2(midpoint = 0, low = brewer.pal(9,'RdYlBu')[9], mid  = brewer.pal(9,'RdYlBu')[5], high = brewer.pal(9,'RdYlBu')[1], na.value = "#FFFFFF", name = as.character(i)) +
    geom_text_repel(data = tmpdata, aes(label = fam_members), size = 4, fontface = "italic", color = "#005500", box.padding = .4, min.segment.length = unit(0.01,'lines'), segment.size = 0.3) +
    xlim(c(-xmax,xmax)) +
    ylim(c(ymin,ymax)) +
    plain.theme})
names(toplot.plots.svn) <- toplot
plot1 <- toplot.plots.svp$Cluster_4_svp + toplot.plots.svp$stable_int_svp + plot_layout(nrow = 2)
plot2 <- toplot.plots.svn$Cluster_4_svn + toplot.plots.svn$stable_int_svn + plot_layout(nrow = 2)
plot3 <- toplot.plots.svp$Cluster_1_svp + toplot.plots.svp$neg_int_svp + plot_layout(nrow = 2)
plot4 <- toplot.plots.svn$Cluster_1_svn + toplot.plots.svn$neg_int_svn + plot_layout(nrow = 2)
#ggsave(plot = plot2, filename = './Pres_plots/increasing_svn.pdf', height = 8, width = 6)
#ggsave(plot = plot3, filename = './Pres_plots/decreasing_svp.pdf', height = 8, width = 6)
currpeaks <- atac.results$peak[which(atac.results$Gene == 'Il10')]
currpeaks <- intersect(currpeaks, rownames(pbm))
temppeaks <- rbind(pbm[currpeaks,])
tfs <- unique(unlist(lapply(1:nrow(temppeaks),function(j){colnames(temppeaks)[which(temppeaks[j,] > 0)]})))
tgttfs <- tfs[which(tfs %in% sig.res$TFs)]
Il10.tfs <- tgttfs
rm(currpeaks,temppeaks,tfs,tgttfs)
x.breaks <- sort(unique(c(-FromZero(values = abs(plotdata$svn_coefs), nbreaks = 2:4),FromZero(values = abs(plotdata$svn_coefs), nbreaks = 2:4))))
y.breaks <- sort(unique(c(-FromZero(values = abs(plotdata$svp_coefs), nbreaks = 2:4),FromZero(values = abs(plotdata$svp_coefs), nbreaks = 2:4))))
plotdata$svn_sig <- plotdata$svn_pvals > 3
plotdata$svp_sig <- plotdata$svp_pvals > 3
plot5 <- ggplot(plotdata, aes(svn_coefs, svp_coefs)) +
  geom_hline(yintercept = 0, linetype = 2, size = 0.5) +
  geom_vline(xintercept = 0, linetype = 2, size = 0.5) +
  geom_point(aes(fill = svn_sig, color = svp_sig), shape = 21, size = 3, stroke = 1) +
  geom_text_repel(data = subset(co1, svn_pvals > 3 | svp_pvals > 3), aes(label = fam_members), size = 4, fontface = "italic", color = "#005500", box.padding = .4, min.segment.length = unit(0.01,'lines'), segment.size = 0.5) +
  scale_x_continuous(limits = c(min(x.breaks),max(x.breaks)), breaks = x.breaks, expand = expansion()) +
  scale_y_continuous(limits = c(min(y.breaks),max(y.breaks)), breaks = y.breaks, expand = expansion()) +
  scale_color_manual(values = c('#8E8E8E','#CD0000')) +
  scale_fill_manual(values = c('#8E8E8E','#000000')) +
  plain.theme +
  theme(legend.position = 'none')
ggsave(plot = plot5, filename = './Pres_plots/coefs_col_all.pdf', width = 5.6, height = 5, dpi = 300)
plot6 <- ggplot(plotdata, aes(svn_coefs, svp_coefs)) +
  geom_hline(yintercept = 0, linetype = 2, size = 0.5) +
  geom_vline(xintercept = 0, linetype = 2, size = 0.5) +
  geom_point(aes(fill = svp_pvals), shape = 21, size = 3) +
  scale_fill_gradientn(colors = brewer.pal(9, 'YlOrRd')) +
  geom_text_repel(data = subset(co1, svp_pvals > 3), aes(label = fam_members), size = 4, fontface = "italic", color = "#005500", box.padding = .4, min.segment.length = unit(0.01,'lines'), segment.size = 0.3) +
  xlim(-0.16,0.16) +
  ylim(-0.00013,0.00013) +
  plain.theme
#ggsave(plot = plot6, filename = './Pres_plots/coefs_col_svp.pdf', width = 6, height = 5, dpi = 300)
hm.data <- read.csv('for_hm.csv', stringsAsFactors = FALSE)
hm.data$Category <- factor(hm.data$Category)
cutoff <- unique(as.vector(na.omit(unlist(lapply(1:nrow(hm.data),function(i){
  if(hm.data$Category[i] == 'TFs'){
    cutoff <- rna.results$Gene[which(rna.results$padj_svn < 0.05)]
  }
  else{
    cutoff <- rna.results$Gene[which(abs(rna.results$lg2FC_svn) > 1 & rna.results$padj_svn < 0.05)]
  }
  cutoff[match(hm.data$Gene[i], cutoff)]
})))))
hm.data <- subset(hm.data, Gene %in% cutoff)
rna.fpkm.zscore <- data.frame(do.call(rbind,lapply(1:nrow(fpkm.counts),function(i){
  as.vector(scale(as.vector(fpkm.counts[i,])))
  })))
colnames(rna.fpkm.zscore) <- colnames(fpkm.counts)
rownames(rna.fpkm.zscore) <- rownames(fpkm.counts)
rna.fpkm.zscore$Gene <- rownames(rna.fpkm.zscore)
hm.data <- merge(hm.data, rna.fpkm.zscore, by.x = 'Gene', by.y = 'Gene', sprt = FALSE)
hm.data <- data.frame(Gene = hm.data$Gene, Category = hm.data$Category, Sample = rep(colnames(fpkm.counts), each = nrow(hm.data)), FPKM = as.vector(as.matrix(hm.data[,grep('RNA',colnames(hm.data))])))
hm.data <- hm.data[-grep('pos', hm.data$Sample),]
hmplots <- lapply(levels(hm.data$Category),function(i){
  curcat <- as.character(i)
  tempdata <- droplevels(subset(hm.data, Category == curcat))
  geneorder <- tempdata[grep('neg',tempdata$Sample),]
  geneorder <- data.frame(Gene = names(tapply(geneorder$FPKM, geneorder$Gene, mean, na.rm = TRUE)), Expr = tapply(geneorder$FPKM, geneorder$Gene, mean, na.rm = TRUE))
  geneorder <- geneorder[order(geneorder$Expr, decreasing = TRUE),]
  geneorder <- as.character(geneorder$Gene)
  ggplot(data = tempdata, aes(Sample,Gene)) +
    geom_tile(aes(fill = FPKM), color = "#FFFFFF", size = 0.2) +
    scale_fill_gradientn(colors = rev(brewer.pal(9,'RdYlBu')), limits = c(min(hm.data$FPKM),max(hm.data$FPKM))) +
    scale_x_discrete(limits = c('RNA_G1_neg','RNA_G2_neg','RNA_G3_neg','RNA_G1_stable','RNA_G2_stable','RNA_G3_stable')) +
    scale_y_discrete(limits = geneorder) +
    plain.theme +
    theme(axis.text.y = element_text(face = 'italic'), axis.text.x = element_blank())
  ggsave(paste0('./Pres_plots/',curcat,'_hm.pdf'), width = 4, height = (0.2*length(geneorder) + 0.5))
})
names(hmplots) <- levels(hm.data$Category)
rna.fpkm.log <- data.frame(do.call(cbind,lapply(1:ncol(fpkm.counts),function(i){log2(fpkm.counts[,i] + 1)})))
colnames(rna.fpkm.log) <- colnames(fpkm.counts)
rownames(rna.fpkm.log) <- rownames(fpkm.counts)
rna.clustersizes$neg_exp_raw <- unlist(lapply(1:nrow(rna.clustersizes),function(i){
  currclust <- rna.clustersizes$Cluster[i]
  if(currclust == 'all'){mean(unlist(rna.fpkm.log[,grep('neg', colnames(rna.fpkm.log))]))}
  else{mean(unlist(rna.fpkm.log[rna.results$Gene[which(rna.results$RNA.cluster == currclust)],grep('neg', colnames(rna.fpkm.log))]))}
}))
rna.clustersizes$pos_exp_raw <- unlist(lapply(1:nrow(rna.clustersizes),function(i){
  currclust <- rna.clustersizes$Cluster[i]
  if(currclust == 'all'){mean(unlist(rna.fpkm.log[,grep('pos', colnames(rna.fpkm.log))]))}
  else{mean(unlist(rna.fpkm.log[rna.results$Gene[which(rna.results$RNA.cluster == currclust)],grep('pos', colnames(rna.fpkm.log))]))}
}))
rna.clustersizes$stable_exp_raw <- unlist(lapply(1:nrow(rna.clustersizes),function(i){
  currclust <- rna.clustersizes$Cluster[i]
  if(currclust == 'all'){mean(unlist(rna.fpkm.log[,grep('stable', colnames(rna.fpkm.log))]))}
  else{mean(unlist(rna.fpkm.log[rna.results$Gene[which(rna.results$RNA.cluster == currclust)],grep('stable', colnames(rna.fpkm.log))]))}
}))
saveRDS(rna.clustersizes,'rna_clustersizes_forHM.rds')
