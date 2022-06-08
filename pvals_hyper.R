library(LimitBreaks)
pvfun <- function(vals){
  unlist(lapply(vals,function(j){
    if(j > 0.05){return('ns')}
    else{
      if(j > 0.01){return('*')}
      else{
        if(j > 0.001){return('**')}
        else{return('***')}
      }
    }
  }))
}
Il10gut$plot_id <- interaction(Il10gut$seurat_clusters, Il10gut$SortID, drop = TRUE, sep = "_")
#elim <- sapply(1:nrow(scrna.df),function(i){
#  length(which(Il10gut$plot_id == scrna.df$scrna_cluster[i]))
#})
#elim <- which(elim > 2)
scrna.df <- scrna.df[which(scrna.df$scrna_cluster %in% c('3_Il10neg','1_Il10neg','2_Il10neg','7_Il10neg','1_Il10pos','2_Il10pos','3_Il10pos','8_Il10neg','8_Il10pos','0_Il10pos','10_Il10pos','9_Il10neg','9_Il10pos')),]
temp1 <- data.frame(clust1 = 1:7, clust2 = c('bI','bII','bIII','bIV','bV','bVI','bVII'))
bulk.clusters$Cluster <- sapply(1:nrow(bulk.clusters),function(i){
  temp1$clust2[which(temp1$clust1 == bulk.clusters$Cluster[i])]
})
scrna.df$P_val_hyper <- unlist(lapply(1:nrow(scrna.df),function(i){
  clust.genes <- rownames(Il10gut@assays$RNA@counts)[which(rowSums(Il10gut@assays$RNA@counts[,which(Il10gut$plot_id == scrna.df$scrna_cluster[i])]) > 0)]
  df0 <- FindMarkers(Il10gut, ident.1 = scrna.df$scrna_cluster[i], group.by = 'plot_id', logfc.threshold = 0.5, min.pct = .1)
  tm <- intersect(bulk.clusters$Gene[which(bulk.clusters$Cluster == scrna.df$bulk_cluster[i])],clust.genes)
  tn <- setdiff(clust.genes,tm)
  tk <- rownames(df0)[which(df0$p_val_adj < 0.05 & df0$avg_log2FC > 0)]
  tq <- intersect(tm,tk)
  phyper(length(tq),length(tm),length(tn),length(tk), lower.tail = FALSE)
}))
scrna.df$P_val_hyper_adj <- p.adjust(scrna.df$P_val_hyper, method = 'holm')
scrna.df$P_val_hyper_adj <- pvfun(scrna.df$P_val_hyper_adj)
scrna.df2 <- data.frame(scrna_cluster = scrna.df$scrna_cluster, do.call(rbind,strsplit(as.character(scrna.df$scrna_cluster), split = '_')))
pv_cols <- data.frame(pval = c('ns','*','**','***'), cols = c("#888888",brewer.pal(9,'BuPu')[2],brewer.pal(9,'BuPu')[4],brewer.pal(9,'BuPu')[6]))
colnames(scrna.df2) <- c('scrna_cluster','cluster','sort_id')
clust_cols2 <- subset(clust_cols2, cluster %in% scrna.df2$cluster)
p.scrna.pval <- ggplot(data = scrna.df, aes(scrna_cluster, bulk_cluster)) +
  geom_tile(aes(fill = P_val_hyper_adj), color = "#FFFFFF", size = 0.3) +
  scale_fill_manual(values = pv_cols$cols[order(as.character(pv_cols$pval))]) +
  new_scale_fill() +
  geom_tile(data = scrna.df2, aes(x = scrna_cluster, y = 'cluster1', fill = cluster), height = 0.5, color = "#FFFFFF", size = 0.3) + 
  scale_fill_manual(values = clust_cols2$cols[order(as.character(clust_cols2$cluster))]) +
  new_scale_fill() +
  geom_tile(data = scrna.df2, aes(x = scrna_cluster, y = 'cluster2', fill = sort_id), height = 0.5, color = "#FFFFFF", size = 0.3) + 
  scale_fill_manual(values = sort_cols) +
  scale_x_discrete(limits = c('3_Il10neg','1_Il10neg','2_Il10neg','7_Il10neg','1_Il10pos','2_Il10pos','3_Il10pos','8_Il10neg','8_Il10pos','0_Il10pos','10_Il10pos','9_Il10neg','9_Il10pos')) +
  scale_y_discrete(limits = c('cluster1','cluster2','bI','bII','bIII','bIV','bV','bVI','bVII')) +
  plain.theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave('pval_bulk_scrna_newgut.pdf', height = 6, width = 10, dpi = 300)
