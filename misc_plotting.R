library(ggplot2)
library(ggrepel)
library(ggridges)
library(patchwork)
library(RColorBrewer)
plain.theme <- theme(axis.line.x = element_line(color = '#000000', size = 0.5, lineend = 'round'),
                     axis.line.y = element_line(color = '#000000', size = 0.5, lineend = 'round'),
                     axis.ticks = element_line(color = 'black', size = 0.5, lineend = 'round'),
                     axis.ticks.length = unit(4,'points'),
                     panel.background = element_rect(fill = 'white'),
                     text = element_text(size = 8, color = 'black'),
                     axis.text = element_text(size = 8, color = 'black'))
rna.fpkmcounts <- readRDS('rna_fpkm.rds')
rna.normcounts <- readRDS('rna_normcounts.rds')
rna.results <- readRDS('rna_results.rds')
rna.sample.info <- readRDS('rna_sampleinfo.rds')
atac.fpkmcounts <- readRDS('atac_fpkm.rds')
atac.normcounts <- readRDS('atac_normcounts.rds')
atac.results <- readRDS('atac_results.rds')
atac.sample.info <- readRDS('atac_sampleinfo.rds')
pbm <- readRDS('./modeling/pbm.rds')
numtfs <- ncol(pbm)
toexclude <- c('Ighv','Ighj','Ighd','Igkv','Igkj','Iglv','Iglj','Trav','Traj','Trbv','Trbj','Trbd','Trgv','Trgj','Trdv','Trdj','Trdd')
toexclude <- unique(c(unlist(lapply(toexclude,function(i){
  rownames(rna.results)[grep(i,rownames(rna.results))]
}))))
toexclude <- unique(c(toexclude,rownames(rna.fpkmcounts)[which(rowMeans(rna.fpkmcounts) <= mean(rna.fpkmcounts['Cd8a',]))]))
toexclude <- unique(c(toexclude,rownames(rna.fpkmcounts)[which(apply(rna.fpkmcounts,1,median) == 0)]))
rna.results <- rna.results[-match(toexclude,rownames(rna.results)),]
rna.normcounts <- rna.normcounts[-match(toexclude,rownames(rna.normcounts)),]
rna.fpkmcounts <- rna.fpkmcounts[-match(toexclude,rownames(rna.fpkmcounts)),]
rm(toexclude)
sig.genes <- unique(c(which(rna.results$padj_pvn < 0.05),which(rna.results$padj_svn < 0.05),which(rna.results$padj_svp < 0.05)))
z.scorecounts <- do.call(rbind,lapply(seq_along(rownames(rna.normcounts)),function(i){
  meanr <- mean(rna.normcounts[i,])
  sdr <- sd(rna.normcounts[i,])
  if(sdr == 0){
    newrow <- rep(0,ncol(rna.normcounts))
  }
  else{
    newrow <- (rna.normcounts[i,] - meanr)/sdr}
}))
rownames(z.scorecounts) <- rownames(rna.normcounts)
colnames(z.scorecounts) <- colnames(rna.normcounts)
rna.forclust <- z.scorecounts[sig.genes,]
rna.clustcol <- hclust(dist(t(rna.forclust)))
rna.kmeans <- readRDS('rna_kmeans.rds')
rna.clustersizes <- table(rna.kmeans$cluster)
rna.sizes.plot <- rna.clustersizes[1]
for(n in 2:length(rna.clustersizes)){rna.sizes.plot <- c(rna.sizes.plot,sum(rna.clustersizes[1:n]))}
rna.clustersizes <- data.frame(rna.clustersizes)
colnames(rna.clustersizes) <- c('Cluster', 'N_genes')
rna.gene.order <- data.frame(Gene = rownames(rna.forclust)[order(rna.kmeans$cluster)], Cluster = rna.kmeans$cluster[order(rna.kmeans$cluster)])
rna.samp.order <- labels(as.dendrogram(rna.clustcol))
rna.kmeans.data <- data.frame(Gene = rownames(rna.forclust), Sample = rep(colnames(rna.forclust), each = nrow(rna.forclust)), Value = as.vector(rna.forclust))
atac.results$RNA.cluster <- unlist(lapply(seq_along(rownames(atac.results)),function(i){
  if(atac.results$Gene[i] %in% rna.gene.order$Gene){paste0('Cluster_', rna.gene.order$Cluster[which(rna.gene.order == atac.results$Gene[i])])}
  else{paste0('No_cluster')}
}))
rna.results$RNA.cluster <- unlist(lapply(seq_along(rownames(rna.results)),function(i){
  if(rna.results$Gene[i] %in% rna.gene.order$Gene){paste0('Cluster_', rna.gene.order$Cluster[which(rna.gene.order == rna.results$Gene[i])])}
  else{paste0('No_cluster')}
}))
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
rna.plt0 <- ggplot(data = subset(rna.results, RNA.cluster == 'No_cluster'), aes(lg2FC_svn,lg2FC_svp)) +
  geom_point(color = 'grey', alpha = 0.4) + 
  geom_point(data = subset(rna.results, RNA.cluster != 'No_cluster'), aes(lg2FC_svn, lg2FC_svp, color = RNA.cluster), alpha = .8) +
  geom_hline(yintercept = 0, linetype = 2, size = 0.5) +
  geom_vline(xintercept = 0, linetype = 2, size = 0.5) +
  geom_text_repel(data = subset(rna.results, Gene == "Xcl1" | Gene == "Il10"), aes(lg2FC_svn,lg2FC_svp, label = Gene), nudge_x = -1, nudge_y = -1, fontface = 'bold', box.padding = 1, min.segment.length = 0.1) +
  scale_color_manual(limits = c("Cluster_1","Cluster_2","Cluster_3","Cluster_4","Cluster_5","Cluster_6","Cluster_7"),
                     labels = c("Cluster I","Cluster II","Cluster III","Cluster IV","Cluster V","Cluster VI","Cluster VII"),
                     values = brewer.pal(7,'Set2')) +
  xlim(min(rna.results$lg2FC_svn),-min(rna.results$lg2FC_svn)) +
  ylim(min(rna.results$lg2FC_svp),-min(rna.results$lg2FC_svp)) +
  plain.theme +
  theme(legend.position = 'none')
atac.plt0 <- ggplot(data = subset(atac.results, RNA.cluster == 'No_cluster'), aes(lg2FC_svn,lg2FC_svp)) +
  geom_point(color = 'grey', alpha = 0.4) + 
  geom_point(data = subset(atac.results, RNA.cluster != 'No_cluster'), aes(lg2FC_svn,lg2FC_svp, color = RNA.cluster), alpha = .8) +
  geom_hline(yintercept = 0, linetype = 2, size = 0.5) +
  geom_vline(xintercept = 0, linetype = 2, size = 0.5) +
  scale_color_manual(limits = c("Cluster_1","Cluster_2","Cluster_3","Cluster_4","Cluster_5","Cluster_6","Cluster_7"),
                     labels = c("Cluster I","Cluster II","Cluster III","Cluster IV","Cluster V","Cluster VI","Cluster VII"),
                     values = brewer.pal(7,'Set2')) +
  xlim(-max(atac.results$lg2FC_svn),max(atac.results$lg2FC_svn)) +
  ylim(min(atac.results$lg2FC_svp),-min(atac.results$lg2FC_svp)) +
  plain.theme +
  theme(legend.position = 'none')
#ggsave(paste0('./Pres_plots/all_clust_rna_atac_fcfc.pdf'), plot = (rna.plt0 + atac.plt0), height = 5, width = 10.5, dpi = 300)
lm.clust.inc <- summary(lm(lg2FC_svp~lg2FC_svn, data = subset(rna.results, RNA.cluster == clust.inc)))
lm.clust.dec <- summary(lm(lg2FC_svp~lg2FC_svn, data = subset(rna.results, RNA.cluster == clust.dec)))
lm.all <- summary(lm(lg2FC_svp~lg2FC_svn, data = rna.results))
temp2 <- data.frame(mod = rep(c('all',clust.inc,clust.dec), each = 2),
                    xval = rep(c(max(abs(rna.results$lg2FC_svn))-6,max(abs(rna.results$lg2FC_svn))), times = 3),
                    lab = c(paste0('R2 = ',format(round(lm.all$adj.r.squared,3), nsmall = 3)),
                            paste0('Slope = ',format(round(lm.all$coefficients[2,1],3), nsmall = 3)),
                            paste0('R2 = ',format(round(lm.clust.inc$adj.r.squared,3), nsmall = 3)),
                            paste0('Slope = ',format(round(lm.clust.inc$coefficients[2,1],3), nsmall = 3)),
                            paste0('R2 = ',format(round(lm.clust.dec$adj.r.squared,3), nsmall = 3)),
                            paste0('Slope = ',format(round(lm.clust.dec$coefficients[2,1],3), nsmall = 3))),
                    yval = rep(c(max(abs(rna.results$lg2FC_svp)),max(abs(rna.results$lg2FC_svp))-1,max(abs(rna.results$lg2FC_svp))-2), each = 2),
                    cols = rep(c('#888888','#d4653a','#3e9a7d'),each = 2))
rna.plt1 <- ggplot(data = subset(rna.results, RNA.cluster != clust.dec & RNA.cluster != clust.inc), aes(lg2FC_svn,lg2FC_svp)) +
  geom_point(color = 'grey', alpha = 0.2) + 
  geom_point(data = subset(rna.results, RNA.cluster == clust.dec | RNA.cluster == clust.inc), aes(lg2FC_svn,lg2FC_svp, color = RNA.cluster), alpha = 0.8) +
  geom_smooth(data = subset(rna.results, RNA.cluster == clust.dec), color = '#3e9a7d', method='lm', formula= y~x) +
  geom_smooth(data = subset(rna.results, RNA.cluster == clust.inc), color = '#d4653a', method='lm', formula= y~x) +
  geom_smooth(data = rna.results, color = 'grey', method='lm', formula= y~x) +
  geom_hline(yintercept = 0, linetype = 2, size = 0.5) +
  geom_vline(xintercept = 0, linetype = 2, size = 0.5) +
  geom_text(data = temp2, aes(x = xval, y = yval, label = lab, color = mod), hjust = 'right') +
  scale_color_manual(limits = c(clust.dec,clust.inc,'all'), labels = c("Cluster I","Cluster IV","all"), values = c('#3e9a7d','#d4653a','#888888')) +
  xlim(-max(temp2$xval),max(temp2$xval)) +
  ylim(-max(temp2$yval),max(temp2$yval)) +
  plain.theme +
  theme(legend.position = 'none')
lm.clust.inc <- summary(lm(lg2FC_svp~lg2FC_svn, data = subset(atac.results, RNA.cluster == clust.inc)))
lm.clust.dec <- summary(lm(lg2FC_svp~lg2FC_svn, data = subset(atac.results, RNA.cluster == clust.dec)))
lm.all <- summary(lm(lg2FC_svp~lg2FC_svn, data = atac.results))
temp2 <- data.frame(mod = rep(c('all',clust.inc,clust.dec), each = 2),
                    xval = rep(c(max(abs(atac.results$lg2FC_svn))-3,max(abs(atac.results$lg2FC_svn))), times = 3),
                    lab = c(paste0('R2 = ',format(round(lm.all$adj.r.squared,3), nsmall = 3)),
                            paste0('Slope = ',format(round(lm.all$coefficients[2,1],3), nsmall = 3)),
                            paste0('R2 = ',format(round(lm.clust.inc$adj.r.squared,3), nsmall = 3)),
                            paste0('Slope = ',format(round(lm.clust.inc$coefficients[2,1],3), nsmall = 3)),
                            paste0('R2 = ',format(round(lm.clust.dec$adj.r.squared,3), nsmall = 3)),
                            paste0('Slope = ',format(round(lm.clust.dec$coefficients[2,1],3), nsmall = 3))),
                    yval = rep(c(max(abs(atac.results$lg2FC_svp)),max(abs(atac.results$lg2FC_svp))-0.25,max(abs(atac.results$lg2FC_svp))-0.5), each = 2),
                    cols = rep(c('#888888','#d4653a','#3e9a7d'),each = 2))
atac.plt1 <- ggplot(data = subset(atac.results, RNA.cluster != clust.dec & RNA.cluster != clust.inc), aes(lg2FC_svn,lg2FC_svp)) +
  geom_point(color = 'grey', alpha = 0.2) + 
  geom_point(data = subset(atac.results, RNA.cluster == clust.dec | RNA.cluster == clust.inc), aes(lg2FC_svn,lg2FC_svp, color = RNA.cluster), alpha = 0.8) +
  geom_smooth(data = subset(atac.results, RNA.cluster == clust.dec), color = '#3e9a7d', method='lm', formula= y~x) +
  geom_smooth(data = subset(atac.results, RNA.cluster == clust.inc), color = '#d4653a', method='lm', formula= y~x) +
  geom_smooth(data = atac.results, color = 'grey', method='lm', formula= y~x) +
  geom_hline(yintercept = 0, linetype = 2, size = 0.5) +
  geom_vline(xintercept = 0, linetype = 2, size = 0.5) +
  geom_text(data = temp2, aes(x = xval, y = yval, label = lab, color = mod), hjust = 'right') +
  scale_color_manual(limits = c(clust.dec,clust.inc,'all'), labels = c("Cluster I","Cluster IV","all"), values = c('#3e9a7d','#d4653a','#888888')) +
  xlim(-max(temp2$xval),max(temp2$xval)) +
  ylim(-max(temp2$yval),max(temp2$yval)) +
  plain.theme +
  theme(legend.position = 'none')
#ggsave(paste0('./Pres_plots/main_clust_rna_atac_fcfc.pdf'), plot = (rna.plt1 + atac.plt1), height = 5, width = 10.5, dpi = 300)
atac.results.v2 <- data.frame(do.call(cbind,list(tapply(atac.results$lg2FC_svn, atac.results$Gene, mean),tapply(atac.results$lg2FC_svp, atac.results$Gene, mean))))
colnames(atac.results.v2) <- c('lg2FC_svn','lg2FC_svp')
atac.results.v2$Gene <- rownames(atac.results.v2)
both.res <- merge(rna.results[,c(1,5,7,11)], atac.results.v2, by = 'Gene', sort = FALSE, suffixes = c('.rna','.atac'))
lm.clust.inc <- summary(lm(lg2FC_svn.atac~lg2FC_svn.rna, data = subset(both.res, RNA.cluster == clust.inc)))
lm.clust.dec <- summary(lm(lg2FC_svn.atac~lg2FC_svn.rna, data = subset(both.res, RNA.cluster == clust.dec)))
lm.all <- summary(lm(lg2FC_svn.atac~lg2FC_svn.rna, data = both.res))
temp2 <- data.frame(mod = rep(c('all',clust.inc,clust.dec), each = 2),
                    xval = rep(c(max(abs(both.res$lg2FC_svn.rna))-5.5,max(abs(both.res$lg2FC_svn.rna))), times = 3),
                    lab = c(paste0('R2 = ',format(round(lm.all$adj.r.squared,3), nsmall = 3)),
                            paste0('Slope = ',format(round(lm.all$coefficients[2,1],3), nsmall = 3)),
                            paste0('R2 = ',format(round(lm.clust.inc$adj.r.squared,3), nsmall = 3)),
                            paste0('Slope = ',format(round(lm.clust.inc$coefficients[2,1],3), nsmall = 3)),
                            paste0('R2 = ',format(round(lm.clust.dec$adj.r.squared,3), nsmall = 3)),
                            paste0('Slope = ',format(round(lm.clust.dec$coefficients[2,1],3), nsmall = 3))),
                    yval = rep(c(max(abs(both.res$lg2FC_svn.atac)),max(abs(both.res$lg2FC_svn.atac))-0.3,max(abs(both.res$lg2FC_svn.atac))-0.6), each = 2),
                    cols = rep(c('#888888','#d4653a','#3e9a7d'),each = 2))
both.plt1 <- ggplot(data = subset(both.res, RNA.cluster != clust.dec & RNA.cluster != clust.inc), aes(lg2FC_svn.rna,lg2FC_svn.atac)) +
  geom_point(color = 'grey', alpha = 0.2) + 
  geom_point(data = subset(both.res, RNA.cluster == clust.dec | RNA.cluster == clust.inc), aes(lg2FC_svn.rna, lg2FC_svn.atac, color = RNA.cluster), alpha = 0.8) +
  geom_smooth(data = subset(both.res, RNA.cluster == clust.dec), color = '#3e9a7d', method='lm', formula= y~x) +
  geom_smooth(data = subset(both.res, RNA.cluster == clust.inc), color = '#d4653a', method='lm', formula= y~x) +
  geom_smooth(data = both.res, color = 'grey', method='lm', formula= y~x) +
  geom_hline(yintercept = 0, linetype = 2, size = 0.5) +
  geom_vline(xintercept = 0, linetype = 2, size = 0.5) +
  geom_text(data = temp2, aes(x = xval, y = yval, label = lab, color = mod), hjust = 'right') +
  scale_color_manual(limits = c(clust.dec,clust.inc,'all'), labels = c("Cluster I","Cluster IV","all"), values = c('#3e9a7d','#d4653a','#888888')) +
  xlim(-max(temp2$xval),max(temp2$xval)) +
  ylim(-max(temp2$yval),max(temp2$yval)) +
  plain.theme +
  theme(legend.position = 'none')
lm.clust.inc <- summary(lm(lg2FC_svp.atac~lg2FC_svp.rna, data = subset(both.res, RNA.cluster == clust.inc)))
lm.clust.dec <- summary(lm(lg2FC_svp.atac~lg2FC_svp.rna, data = subset(both.res, RNA.cluster == clust.dec)))
lm.all <- summary(lm(lg2FC_svp.atac~lg2FC_svp.rna, data = both.res))
temp2 <- data.frame(mod = rep(c('all',clust.inc,clust.dec), each = 2),
                    xval = rep(c(max(abs(both.res$lg2FC_svp.rna))-5,max(abs(both.res$lg2FC_svp.rna))), times = 3),
                    lab = c(paste0('R2 = ',format(round(lm.all$adj.r.squared,3), nsmall = 3)),
                            paste0('Slope = ',format(round(lm.all$coefficients[2,1],3), nsmall = 3)),
                            paste0('R2 = ',format(round(lm.clust.inc$adj.r.squared,3), nsmall = 3)),
                            paste0('Slope = ',format(round(lm.clust.inc$coefficients[2,1],3), nsmall = 3)),
                            paste0('R2 = ',format(round(lm.clust.dec$adj.r.squared,3), nsmall = 3)),
                            paste0('Slope = ',format(round(lm.clust.dec$coefficients[2,1],3), nsmall = 3))),
                    yval = rep(c(max(abs(both.res$lg2FC_svp.atac)),max(abs(both.res$lg2FC_svp.atac))-0.1,max(abs(both.res$lg2FC_svp.atac))-0.2), each = 2),
                    cols = rep(c('#888888','#d4653a','#3e9a7d'),each = 2))
both.plt2 <- ggplot(data = subset(both.res, RNA.cluster != clust.dec & RNA.cluster != clust.inc), aes(lg2FC_svp.rna,lg2FC_svp.atac)) +
  geom_point(color = 'grey', alpha = 0.2) + 
  geom_point(data = subset(both.res, RNA.cluster == clust.dec | RNA.cluster == clust.inc), aes(lg2FC_svp.rna, lg2FC_svp.atac, color = RNA.cluster), alpha = 0.8) +
  geom_smooth(data = subset(both.res, RNA.cluster == clust.dec), color = '#3e9a7d', method='lm', formula= y~x) +
  geom_smooth(data = subset(both.res, RNA.cluster == clust.inc), color = '#d4653a', method='lm', formula= y~x) +
  geom_smooth(data = both.res, color = 'grey', method='lm', formula= y~x) +
  geom_hline(yintercept = 0, linetype = 2, size = 0.5) +
  geom_vline(xintercept = 0, linetype = 2, size = 0.5) +
  geom_text(data = temp2, aes(x = xval, y = yval, label = lab, color = mod), hjust = 'right') +
  scale_color_manual(limits = c(clust.dec,clust.inc,'all'), labels = c("Cluster I","Cluster IV","all"), values = c('#3e9a7d','#d4653a','#888888')) +
  xlim(-max(temp2$xval),max(temp2$xval)) +
  ylim(-max(temp2$yval),max(temp2$yval)) +
  plain.theme +
  theme(legend.position = 'none')
#ggsave(plot = (both.plt1 + both.plt2), filename = paste0('./Pres_plots/rna_atac_','int_clusts','.pdf'), width = 10.5, height = 5)
atac.density <- do.call(cbind,lapply(1:ncol(atac.fpkmcounts),function(i){log2(atac.fpkmcounts[,i] + 1)}))
colnames(atac.density) <- colnames(atac.normcounts)
#atac.density <- atac.normcounts
testaves <- do.call(rbind,lapply(seq_along(rownames(atac.density)),function(i){
  tapply(as.vector(atac.density[i,]),atac.sample.info$Celltype,mean)
}))
colnames(testaves) <- c('neg','pos','stable')
rownames(testaves) <- rownames(atac.density)
atac.density <- data.frame(Celltype = rep(colnames(testaves), each = nrow(testaves)), peak = rep(rownames(testaves), times = 3), stringsAsFactors = FALSE)
atac.density$Value <- as.vector(testaves)
atac.density$Cluster <- atac.results$RNA.cluster[match(atac.density$peak,atac.results$peak)]
plot1 <- ggplot(data = subset(atac.density, Cluster == clust.dec | Cluster == clust.inc), aes(Value,Celltype)) +
  facet_grid(Cluster~.) +
  geom_density_ridges(aes(fill = Celltype), alpha = 1) +
  scale_fill_manual(values = c("#56b4e9", "#d55e00","#e69f00")) +
  #geom_vline(xintercept = density(atac.density$Value[which(atac.density$Celltype == 'stable')])$x[which.max(density(atac.density$Value[which(atac.density$Celltype == 'stable')])$y)], color = "#ff8200") +
  geom_vline(xintercept = testaves['peak_2961','stable']) +
  geom_vline(xintercept = testaves['peak_2961','neg']) +
  scale_x_continuous(limits = c(0.5,7.8), expand = expansion(mult = c(0,0), add = c(0,0))) + 
  plain.theme
#ggsave('./Pres_plots/atac_clust_densities.pdf', height = 6, width = 4, dpi = 300)

tcr.comp <- c('Cd247','Cd4','Cd3e','Tcrb')
tcr.pos <- c('Tespa1','Themis','Plcg1','Pik3r5','Map3k14','Map3k8','Prkcq','Tec','Malt1','Nck2','Itk','Txk')
tcr.neg <- c('Ubash3a','Ubash3b','Dusp6','Cd5','Lag3','Cradd')
tcr.list <- list(tcr.comp,tcr.pos,tcr.neg)
names(tcr.list) <- c('tcr_comp','tcr_pos','tcr_neg')
rna.fpkm.log <- do.call(cbind,lapply(1:ncol(rna.fpkmcounts),function(i){log2(rna.fpkmcounts[,i] + 1)}))
tcrplots <- lapply(tcr.list,function(i){
  geneorder <- rna.results[i,][order(rna.results[i,]$lg2FC_svn),]$Gene
  plotdf <- data.frame(counts = as.vector(as.matrix(t(rna.fpkm.log[i,]))),
             Celltype = rna.sample.info$Celltype,
             gene = rep(i, each = ncol(rna.fpkm.log)))
  plotlims <- testfun1(plotdf$counts)
  pvaldf <- data.frame(gene = rna.results[i,]$Gene, svn = round(rna.results[i,]$padj_svn, 5), svp = round(rna.results[i,]$padj_svp, 5), pvn = round(rna.results[i,]$padj_pvn, 5))
  pvaldf$yval <- unlist(lapply(pvaldf$gene,function(j){
    min(max(plotdf$counts[which(plotdf$gene == j)])*1.1,max(plotlims)*.95)
  }))
  currentplot <- ggplot(data = plotdf, aes(gene, counts)) +
    stat_summary(aes(fill = Celltype), size = 0.5, fun = 'mean', geom = 'col', width = 0.7, position = position_dodge2(), color = "#000000") +
    geom_point(aes(fill = Celltype), shape = 21, size = 2, position = position_jitterdodge(dodge.width = 0.7, jitter.width = 0.2)) +
    geom_text(data = pvaldf, aes(gene, yval, label = svn), color = "#e69f00", nudge_x = -0.2, size = 2) +
    geom_text(data = pvaldf, aes(gene, yval, label = svp), color = "#d55e00", nudge_x = 0.2, size = 2) +
    geom_text(data = pvaldf, aes(gene, yval, label = pvn), color = "#56b4e9", nudge_x = 0, size = 2) +
    scale_fill_manual(values = c("#56b4e9", "#d55e00","#e69f00")) +
    scale_x_discrete(limits = geneorder) +
    scale_y_continuous(limits = c(min(plotlims),max(plotlims)), breaks = plotlims, expand = expansion(mult = c(0,0), add = c(0,0))) +
    plain.theme +
    theme(legend.position = 'none')
})
ggsave('./Pres_plots/tcr_comp.pdf', plot = tcrplots$tcr_comp, height = 2.5, width = (0.5+(length(tcr.list$tcr_comp)*.6)), dpi = 300)
ggsave('./Pres_plots/tcr_pos.pdf', plot = tcrplots$tcr_pos, height = 2.5, width = (0.5+(length(tcr.list$tcr_pos)*.6)), dpi = 300)
ggsave('./Pres_plots/tcr_neg.pdf', plot = tcrplots$tcr_neg, height = 2.5, width = (0.5+(length(tcr.list$tcr_neg)*.6)), dpi = 300)

genesets <- list.files(path = "./Genesets/", pattern = '.txt')
genesets <- unlist(lapply(seq_along(genesets),function(i){
  cgs <- genesets[i]
  substring(cgs,1,regexpr('.txt',cgs)-1)
}))
genesetnames <- unlist(lapply(strsplit(genesets, split = "\\."),paste, collapse = "_"))
genesets <- lapply(genesets,function(i){
  tmp <- read.delim(paste0('./Genesets/',as.character(i),'.txt'), header = TRUE, stringsAsFactors = FALSE, sep = " ")[,2]
})
names(genesets) <- genesetnames
dir.create(paste0('./Pres_plots/CDFs/'), showWarnings = FALSE)
comp.combn <- data.frame(combn(rev(unique(as.character(rna.sample.info$Celltype))),2))
CDF.plots2 <- lapply(seq_along(genesets),function(i){
  cfam <- names(genesets)[i]
  famgenes <- intersect(unique(genesets[[i]]),rownames(rna.results))
  cdfplots <- lapply(comp.combn,function(j){
    temp.a <- paste0(substring(as.character(j[1]),1,1),'v',substring(as.character(j[2]),1,1))
    dir.create(paste0('./Pres_plots/CDFs/',temp.a,'/'), showWarnings = FALSE)
    temp.b <- intersect(grep('lg2FC', colnames(rna.results)),grep(temp.a, colnames(rna.results)))
    data.cdf <- data.frame(FC = c(rna.results[,temp.b],
                                  rna.results[famgenes,temp.b]),
                           Group = factor(c(rep('All', nrow(rna.results)),
                                            rep(cfam, length(famgenes))),
                                          levels = c('All',cfam)))
    cdf.plot <- ggplot(data = data.cdf, aes(x = FC, color = Group)) +
      stat_ecdf(geom = 'step', size = .5) +
      labs(y = 'Cumulative fraction', x = paste0('log2 FC ', temp.a)) +
      xlim(-4,4) +
      scale_color_manual(values = c('#000000','#AA0000'), labels = c(paste0('All genes (n=', nrow(rna.results),')'),
                                                                     paste0(cfam,' (n= ',length(famgenes),')\nP-val: ', formatC(t.test(rna.results[famgenes,temp.b],rna.results[,temp.b])$p.value, format= 'e', digits = 3),'\nP-val: ', formatC(ks.test(rna.results[famgenes,temp.b],rna.results[,temp.b])$p.value, format= 'e', digits = 3))), guide = guide_legend(title = NULL)) +
      plain.theme +
      theme(plot.margin = unit(c(0,1.3,0,0), 'in'), legend.position = c(1.1,.4))
    ggsave(paste0('./Pres_plots/CDFs/',temp.a,'/RNA_',cfam,'.pdf'), height = 2, width = 4.5, dpi = 300)
  })
})
CDF.plots3 <- lapply(c('moore','agl','both'),function(i){
  cfam <- as.character(i)
  famgenes.up <- intersect(grep(cfam, names(genesets)),grep('up', names(genesets)))
  famgenes.up <- intersect(unique(genesets[[famgenes.up]]),rownames(rna.results))
  famgenes.down <- intersect(grep(cfam, names(genesets)),grep('down', names(genesets)))
  famgenes.down <- intersect(unique(genesets[[famgenes.down]]),rownames(rna.results))
  cdfplots <- lapply(comp.combn,function(j){
    temp.a <- paste0(substring(as.character(j[1]),1,1),'v',substring(as.character(j[2]),1,1))
    dir.create(paste0('./Pres_plots/CDFs/',temp.a,'/'), showWarnings = FALSE)
    temp.b <- intersect(grep('lg2FC', colnames(rna.results)),grep(temp.a, colnames(rna.results)))
    data.cdf <- data.frame(FC = c(rna.results[,temp.b],
                                  rna.results[famgenes.up,temp.b],
                                  rna.results[famgenes.down,temp.b]),
                           Group = factor(c(rep('All', nrow(rna.results)),
                                            rep('up', length(famgenes.up)),
                                            rep('down', length(famgenes.down))),
                                          levels = c('All','up','down')))
    cdf.plot <- ggplot(data = data.cdf, aes(x = FC, color = Group)) +
      stat_ecdf(geom = 'step', size = .5) +
      labs(y = 'Cumulative fraction', x = paste0('log2 FC ', temp.a)) +
      xlim(-4,4) +
      scale_color_manual(values = c('#000000','#770000',"#000077"), labels = c(paste0('All genes (n=', nrow(rna.results),')'),
                                                                     paste0('up',' (n= ',length(famgenes.up),')\nt-test: ', formatC(t.test(rna.results[famgenes.up,temp.b],rna.results[,temp.b])$p.value, format= 'e', digits = 3),'\nks-test: ', formatC(ks.test(rna.results[famgenes.up,temp.b],rna.results[,temp.b])$p.value, format= 'e', digits = 3)),
                                                                     paste0('down',' (n= ',length(famgenes.down),')\nt-test: ', formatC(t.test(rna.results[famgenes.down,temp.b],rna.results[,temp.b])$p.value, format= 'e', digits = 3),'\nks-test: ', formatC(ks.test(rna.results[famgenes.down,temp.b],rna.results[,temp.b])$p.value, format= 'e', digits = 3))), guide = guide_legend(title = NULL)) +
      plain.theme +
      theme(plot.margin = unit(c(0.01,1.1,0.01,0.01), 'in'), legend.position = c(1.1,0.5))
    ggsave(paste0('./Pres_plots/CDFs/',temp.a,'/RNA_',cfam,'.pdf'), height = 2, width = 4.5, dpi = 300)
  })
})
testplots <- lapply(unique(both.res$RNA.cluster),function(i){
  clust.cur <- as.character(i)
  both.pltt1 <- ggplot(data = subset(both.res, RNA.cluster != clust.cur), aes(lg2FC_svn.rna,lg2FC_svn.atac)) +
    geom_point(color = 'grey', alpha = 0.2) + 
    geom_point(data = subset(both.res, RNA.cluster == clust.cur), aes(lg2FC_svn.rna, lg2FC_svn.atac), color = '#ea4311', alpha = 0.8) +
    geom_smooth(data = subset(both.res, RNA.cluster == clust.cur), color = '#ea4311', method='lm', formula= y~x) +
    geom_smooth(data = both.res, color = 'grey', method='lm', formula= y~x) +
    geom_hline(yintercept = 0, linetype = 2, size = 0.5) +
    geom_vline(xintercept = 0, linetype = 2, size = 0.5) +
    xlim(-temp1$xval,temp1$xval) +
    ylim(-temp1$yval,temp1$yval) +
    geom_text(data = temp1, aes(x = xval, y = yval), label = round(cor(subset(both.res, RNA.cluster == clust.cur)$lg2FC_svn.rna, subset(both.res, RNA.cluster == clust.cur)$lg2FC_svn.atac),3), color = '#ea4311', nudge_x = -1, nudge_y = -0.5) +
    geom_text(data = temp1, aes(x = xval, y = yval), label = round(cor(both.res$lg2FC_svn.rna, both.res$lg2FC_svn.atac),3), color = '#888888', nudge_x = -1) +  
    plain.theme +
    theme(legend.position = 'none')
  both.pltt2 <- ggplot(data = subset(both.res, RNA.cluster != clust.cur), aes(lg2FC_svp.rna,lg2FC_svp.atac)) +
    geom_point(color = 'grey', alpha = 0.2) + 
    geom_point(data = subset(both.res, RNA.cluster == clust.cur), aes(lg2FC_svp.rna, lg2FC_svp.atac), color = '#ea4311', alpha = 0.8) +
    geom_smooth(data = subset(both.res, RNA.cluster == clust.cur), color = '#ea4311', method='lm', formula= y~x) +
    geom_smooth(data = both.res, color = 'grey', method='lm', formula= y~x) +
    geom_hline(yintercept = 0, linetype = 2, size = 0.5) +
    geom_vline(xintercept = 0, linetype = 2, size = 0.5) +
    xlim(-temp1$xval,temp1$xval) +
    ylim(-temp1$yval,temp1$yval) +
    geom_text(data = temp1, aes(x = xval, y = yval), label = round(cor(subset(both.res, RNA.cluster == clust.cur)$lg2FC_svp.rna, subset(both.res, RNA.cluster == clust.cur)$lg2FC_svp.atac),3), color = '#ea4311', nudge_x = -1, nudge_y = -0.5) +
    geom_text(data = temp1, aes(x = xval, y = yval), label = round(cor(both.res$lg2FC_svp.rna, both.res$lg2FC_svp.atac),3), color = '#888888', nudge_x = -1) +  
    plain.theme +
    theme(legend.position = 'none')
  ggsave(plot = (both.pltt1 + both.pltt2), filename = paste0('./Pres_plots/rna_atac_',clust.cur,'.pdf'), width = 10.5, height = 5)
})

goi.s <- intersect(rownames(rna.fpkm.log),unique(readRDS('../../mm39/mm_cc_genes.rds')$s.genes))
goi.s <- setdiff(goi.s,'Exo1')
s.genes <- data.frame(counts = as.vector(as.matrix(t(rna.fpkm.log[goi.s,]))),
                      zscore = as.vector(as.matrix(t(z.scorecounts[goi.s,]))),
                      fc_norm = as.vector(t(rna.normcounts[goi.s,]/rna.results[goi.s,'neg_lg2expr'])),
                     Celltype = rna.sample.info$Celltype,
                     Sample = rna.sample.info$Sample,
                     gene = rep(goi.s, each = ncol(rna.fpkm.log)))
p1 <- ggplot(data = s.genes, aes(Sample, counts)) +
  geom_point(aes(fill = Celltype), shape = 21, position = position_jitter(width = 0.4), size = 3) +
  stat_summary(aes(group = Sample), fun = 'mean', geom = 'crossbar', size = 0.2, width = 0.7) +
  scale_fill_manual(values = c("#56b4e9", "#d55e00","#e69f00")) +
  scale_y_continuous(limits = c(0,8), breaks = function(a){seq(a[1],a[2],2)}, expand = expansion(0,0)) +
  scale_x_discrete(limits = c('RNA_G1_neg','RNA_G2_neg','RNA_G3_neg','RNA_G1_pos','RNA_G2_pos','RNA_G3_pos','RNA_G1_stable','RNA_G2_stable','RNA_G3_stable')) +
  plain.theme + 
  theme(legend.position = 'none', plot.title = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(angle = 0, hjust = 0.5))
p2 <- ggplot(data = s.genes, aes(Sample,gene)) +
  geom_tile(aes(fill = zscore), color = "#FFFFFF", size = 0.2) +
  scale_fill_gradientn(colors = rev(brewer.pal(11,'RdYlBu')), limits = c(-3,3)) + 
  scale_y_discrete(limits = unique(subset(s.genes, Celltype == 'pos')$gene[order(subset(s.genes, Celltype == 'pos')$counts)])) + 
  scale_x_discrete(limits = c('RNA_G1_neg','RNA_G2_neg','RNA_G3_neg','RNA_G1_pos','RNA_G2_pos','RNA_G3_pos','RNA_G1_stable','RNA_G2_stable','RNA_G3_stable')) +
  plain.theme + 
  theme(plot.title = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(angle = 0, hjust = 0.5))
p3 <- ggplot(data = s.genes, aes(Sample,gene)) +
  geom_tile(aes(fill = fc_norm), color = "#FFFFFF", size = 0.2) +
  scale_fill_gradientn(colors = rev(hcl.colors(256,'RdBu')), limits = c(0.5,1.5)) + 
  scale_y_discrete(limits = unique(subset(s.genes, Celltype == 'pos')$gene[order(subset(s.genes, Celltype == 'pos')$counts)])) + 
  scale_x_discrete(limits = c('RNA_G1_neg','RNA_G2_neg','RNA_G3_neg','RNA_G1_pos','RNA_G2_pos','RNA_G3_pos','RNA_G1_stable','RNA_G2_stable','RNA_G3_stable')) +
  plain.theme + 
  theme(plot.title = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(angle = 0, hjust = 0.5))
ggsave(plot = p2, filename = './Pres_plots/s_genes_hm.pdf', height = 5, width = 4, dpi = 300)

