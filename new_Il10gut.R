library(dplyr)
library(Seurat)
library(RColorBrewer)
library(ggplot2)
library(ggnewscale)
library(patchwork)
library(Matrix)
library(justFDG)
library(ggridges)
plain.theme <- theme(axis.line.x = element_line(color = '#000000', size = 0.5, lineend = 'round'),
                     axis.line.y = element_line(color = '#000000', size = 0.5, lineend = 'round'),
                     axis.ticks = element_line(color = 'black', size = 0.5, lineend = 'round'),
                     axis.ticks.length = unit(4,'points'),
                     panel.background = element_rect(fill = 'white'),
                     text = element_text(size = 10, color = 'black'),
                     axis.text = element_text(size = 10, color = 'black'))
### can be skipped
samples <- c('E2_Il10neg','E2_Il10pos')
data.ls <- dir('../Rawish_data/')
samples.list <- lapply(samples,function(i){
  tempdata <- Read10X(data.dir = paste0('../Rawish_data/',as.character(i),'/filtered_feature_bc_matrix/'), strip.suffix = TRUE)
  temp <- data.ls[grep(as.character(i),data.ls)]
  if(length(grep('hto',temp)) > 0){
    htodata <- Read10X(data.dir = paste0('../Rawish_data/',temp[grep('hto',temp)],'/umi_count/'), gene.column = 1)
    joint.bcs <- intersect(colnames(tempdata), colnames(htodata))
    tempdata <- tempdata[,joint.bcs]
    htodata <- htodata[,joint.bcs]
    rownames(htodata) <- c('LILP','mesLN','Lung','medLN','unmap')
    tempdata <- CreateSeuratObject(counts = tempdata, project = as.character(i))
    tempdata[["HTO"]] <- CreateAssayObject(counts = htodata)
    tempdata <- NormalizeData(tempdata, assay = "HTO", normalization.method = "CLR")
    tempdata <- HTODemux(tempdata, assay = "HTO", positive.quantile = 0.99999999999, kfunc = 'kmeans')
    tempdata <- subset(tempdata, hash.ID != 'Doublet')
    tempdata <- subset(tempdata, hash.ID != 'Negative')
    tempdata <- subset(tempdata, hash.ID != 'unmap')
  }
  else{
    tempdata <- CreateSeuratObject(counts = tempdata, project = as.character(i))
  }
  tempdata <- NormalizeData(tempdata)
  tempdata <- FindVariableFeatures(tempdata, selection.method = "vst", nfeatures = 2000)
  tempdata[["percent.mt"]] <- PercentageFeatureSet(tempdata, pattern = "^mt-")
  Idents(tempdata) <- 'orig.ident'
  tempdata <- subset(tempdata, percent.mt < 10)
  min.co <- unname(quantile(tempdata$nFeature_RNA, probs = seq(0,1,.02))[2])
  max.co <- unname(quantile(tempdata$nFeature_RNA, probs = seq(0,1,.02))[50])
  tempdata <- subset(tempdata, nFeature_RNA > min.co & nFeature_RNA < max.co)
})
names(samples.list) <- samples
Il10combined <- merge(samples.list$E2_Il10neg, y = samples.list$E2_Il10pos, add.cell.ids = c("Il10neg", "Il10pos"), project = "Il10_E2")
Il10combined$SortID <- substring(Il10combined$orig.ident,4,1000)
Il10combined <- FindVariableFeatures(Il10combined, selection.method = "vst", nfeatures = 2000)
Il10combined <- ScaleData(Il10combined, features = VariableFeatures(object = Il10combined))
Il10combined <- RunPCA(Il10combined, features = VariableFeatures(object = Il10combined))
Il10combined <- FindNeighbors(Il10combined, reduction = 'pca', dims = 1:30)
Il10combined <- FindClusters(Il10combined, resolution = 0.5)
og.markers <- FindAllMarkers(Il10combined, logfc.threshold = 1)
og.markers <- og.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
Il10combined <- subset(Il10combined, seurat_clusters != unique(as.character(og.markers$cluster[grep('Ifit',og.markers$gene)])))
Il10gut <- droplevels(subset(Il10combined, hash.ID %in% c('LILP','mesLN')))
rm(samples.list,Il10combined)
Il10gut <- ScaleData(Il10gut, features = rownames(Il10gut))
Il10gut <- FindVariableFeatures(Il10gut, selection.method = "vst", nfeatures = 3000)
Il10gut <- RunPCA(Il10gut, features = VariableFeatures(object = Il10gut))
Il10gut <- FindNeighbors(Il10gut, reduction = 'pca', dims = 1:30)
Il10gut <- FindClusters(Il10gut, resolution = 0.5)
Il10gut <- RunUMAP(Il10gut, reduction = 'pca', dims = 1:30)
df1 <- runFDG(Il10gut@reductions$pca@cell.embeddings[,1:30],Il10gut@graphs$RNA_snn, python.addr = "python")
colnames(df1) <- c('FDG_1','FDG_2')
Il10gut[['FDG']] <- CreateDimReducObject(as.matrix(df1))
rm(df1)
og.markers <- FindAllMarkers(Il10gut, logfc.threshold = 1)
og.markers <- og.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
t1 <- intersect(rownames(as.data.frame.matrix(as.matrix(table(Il10gut$seurat_clusters,Il10gut$hash.ID))/rowSums(as.matrix(table(Il10gut$seurat_clusters,Il10gut$hash.ID)))))[which(as.data.frame.matrix(as.matrix(table(Il10gut$seurat_clusters,Il10gut$hash.ID))/rowSums(as.matrix(table(Il10gut$seurat_clusters,Il10gut$hash.ID))))$mesLN > .9)],
          unique(as.character(og.markers$cluster[grep('Tcf7',og.markers$gene)])))
t1 <- Cells(subset(Il10gut, seurat_clusters == t1))
Il10gut <- RunHP(object = Il10gut, early.cell = sample(t1,1), n.iter = 600, nPCs = 30)
saveRDS(Il10gut,'Il10newgut_ms.rds')
###
Il10gut <- readRDS('Il10newgut_ms.rds')
og.markers <- FindAllMarkers(Il10gut, logfc.threshold = 1)
og.markers <- og.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
cc.genes.mm <- readRDS('../../mm39/mm_cc_genes.rds')
Il10gut <- CellCycleScoring(Il10gut, s.features = cc.genes.mm$s.genes, g2m.features = cc.genes.mm$g2m.genes, set.ident = TRUE)
clust_cols <- ceiling(length(unique(Il10gut$seurat_clusters))/2)
clust_cols <- cbind(hcl.colors(clust_cols,'Pastel 1'),hcl.colors(clust_cols,'Dark 2'))
clust_cols <- as.vector(t(clust_cols))
clust_cols <- clust_cols[1:length(unique(Il10gut$seurat_clusters))]
sort_cols <- c("#c8c8c8", "#f1623d")
re.lab <- read.csv('re_label.csv')
Il10gut$clust_sort <- droplevels(interaction(Il10gut$seurat_clusters,Il10gut$SortID, sep = "_"))
re.lab$clust_sort <- paste(re.lab$sc_clust,re.lab$sort_id, sep = "_")
Il10gut$bulk_id <- re.lab$bulk_id[match(Il10gut$clust_sort,re.lab$clust_sort)]
bulk.cols <- c("#56b4e9","#d55e00","#e69f00","#CCCCCC","#AAAAAA","#323232")
df2 <- as.data.frame(table(Il10gut$seurat_clusters,Il10gut$SortID))
p0 <- ggplot(data = df2, aes(Var1,Freq)) +
  geom_col(aes(fill = Var2)) +
  scale_fill_manual(values = sort_cols) +
  scale_y_continuous(limits = c(0,6000), expand = expansion()) +
  plain.theme
ggsave('counts_by_cluster.pdf', height = 3, width = 5, dpi = 300)
p3 <- DimPlot(Il10gut, reduction = 'harmony', group.by = 'seurat_clusters', pt.size = 0.3, cols = clust_cols) + theme(axis.text = element_blank(), axis.ticks = element_blank())
ggsave('fdg_clusters_ms.pdf', height = 5, width = 5, dpi = 300)
p2 <- DimPlot(Il10gut, reduction = 'harmony', pt.size = 0.3, cols = sort_cols, group.by = 'SortID') + ylim(ggplot_build(p3)$layout$panel_scales_y[[1]]$range$range) + xlim(ggplot_build(p3)$layout$panel_scales_x[[1]]$range$range) + theme(axis.text = element_blank(), axis.ticks = element_blank())
ggsave('fdg_sort_ms.pdf', height = 5, width = 5, dpi = 300)
p4 <- DimPlot(Il10gut, reduction = 'harmony', pt.size = 0.3, cols = bulk.cols, group.by = 'bulk_id')+ ylim(ggplot_build(p3)$layout$panel_scales_y[[1]]$range$range) + xlim(ggplot_build(p3)$layout$panel_scales_x[[1]]$range$range) + theme(axis.text = element_blank(), axis.ticks = element_blank())
ggsave('fdg_bulk_ms.pdf', height = 5, width = 5, dpi = 300)
p5 <- FeaturePlot(Il10gut, reduction = 'harmony', 'S.Score', pt.size = 0.3, order = TRUE) + scale_color_gradientn(colors = brewer.pal(9,'RdPu')) + ylim(ggplot_build(p3)$layout$panel_scales_y[[1]]$range$range) + xlim(ggplot_build(p3)$layout$panel_scales_x[[1]]$range$range) + theme(axis.text = element_blank(), axis.ticks = element_blank())
ggsave('fdg_s_ms.pdf', height = 5, width = 5, dpi = 300)
p6 <- FeaturePlot(Il10gut, reduction = 'harmony', 'entropy', pt.size = 0.3, order = TRUE) + scale_color_gradientn(colors = brewer.pal(9,'YlGnBu')) + ylim(ggplot_build(p3)$layout$panel_scales_y[[1]]$range$range) + xlim(ggplot_build(p3)$layout$panel_scales_x[[1]]$range$range) + theme(axis.text = element_blank(), axis.ticks = element_blank())
ggsave('fdg_entropy_ms.pdf', height = 5, width = 5, dpi = 300)
p7 <- FeaturePlot(Il10gut, reduction = 'harmony', 'pseudo', pt.size = 0.3, order = TRUE) + scale_color_gradientn(colors = brewer.pal(9,'YlGnBu')) + ylim(ggplot_build(p3)$layout$panel_scales_y[[1]]$range$range) + xlim(ggplot_build(p3)$layout$panel_scales_x[[1]]$range$range) + theme(axis.text = element_blank(), axis.ticks = element_blank())
ggsave('fdg_pseudo.pdf', height = 5, width = 5, dpi = 300)
#p1 + p2 + p3 + p4
pall <- p3 + p2 + p4 + p5 + p6 + plot_layout(ncol = 3)
ggsave('fdg_all_ms.pdf', height = 10, width = 16, dpi = 300)
bulk.clusters <- readRDS('../../Il10_LI_genomics/Final_version/rna_kmeans.rds')
bulk.clusters <- data.frame(bulk.clusters$cluster)
bulk.clusters$Gene <- rownames(bulk.clusters)
bulk.clusters$Cluster <- factor(bulk.clusters$bulk.clusters.cluster)
bulk.clusters$bulk.clusters.cluster <- NULL
bulk.clusters <- bulk.clusters[which(bulk.clusters$Gene %in% rownames(Il10gut)),]
bulk.clusters.list <- lapply(levels(bulk.clusters$Cluster),function(i){intersect(bulk.clusters$Gene[which(bulk.clusters$Cluster == as.character(i))], rownames(Il10gut))})
Il10gut <- AddModuleScore(Il10gut, features = bulk.clusters.list, name = 'Cluster_')
p.vlns <- VlnPlot(Il10gut, features = paste0('Cluster_',1:7), cols = clust_cols, group.by = 'seurat_clusters', pt.size = 0) + theme(plot.title = element_text(size = 5)) + plot_layout(ncol = 3)
ggsave('cluster_violins_forMS.pdf', width = 8.2, height = 8, dpi = 300)
scrna.df <- aggregate(FetchData(Il10gut,c(paste('Cluster',sort(unique(bulk.clusters$Cluster)), sep = "_"),'entropy','pseudo')),
                      by = list(interaction(Il10gut$seurat_clusters, Il10gut$SortID, drop = TRUE, sep = "_")), mean, na.rm = TRUE)
colnames(scrna.df) <- c('Cluster','bI','bII','bIII','bIV','bV','bVI','bVII','Entropy','Pseudo')
scrna.df[,c('bI','bII','bIII','bIV','bV','bVI','bVII')] <- scale(scrna.df[,c('bI','bII','bIII','bIV','bV','bVI','bVII')])
scrna.df <- data.frame(scrna_cluster = rep(scrna.df$Cluster,(ncol(scrna.df)-1)), bulk_cluster = rep(colnames(scrna.df)[2:ncol(scrna.df)], each = nrow(scrna.df)), expr = unlist(scrna.df[,2:ncol(scrna.df)]))
scrna.df2 <- data.frame(scrna_cluster = scrna.df$scrna_cluster, do.call(rbind,strsplit(as.character(scrna.df$scrna_cluster), split = '_')))
clust_cols2 <- data.frame(cols = clust_cols, cluster = 0:10)
colnames(scrna.df2) <- c('scrna_cluster','cluster','sort_id')
rna.fpkm <- readRDS('../../Il10_LI_genomics/Final_version/rna_fpkm.rds')
rna.fpkm <- log2(rna.fpkm + 1)
rna.fpkm <- rna.fpkm[bulk.clusters$Gene,]
rna.df <- do.call(cbind,lapply(1:ncol(rna.fpkm),function(i){
  tapply(rna.fpkm[,i],bulk.clusters$Cluster, mean, na.rm = TRUE)
}))
rna.df <- data.frame(Sample = colnames(rna.fpkm),
                     scale(t(rna.df)))
colnames(rna.df) <- c('Sample','bI','bII','bIII','bIV','bV','bVI','bVII')
rna.df <- data.frame(Celltype = rep(rna.df$Sample,7), bulk_cluster = rep(colnames(rna.df)[2:8], each = nrow(rna.df)), expr = unlist(rna.df[,2:8]))
rna.df2 <- data.frame(Celltype = rna.df$Celltype, sort_pop = do.call(rbind,strsplit(as.character(rna.df$Celltype), split = '_'))[,3])
p.scrna <- ggplot(data = scrna.df, aes(scrna_cluster, bulk_cluster)) +
  geom_tile(aes(fill = expr), color = "#FFFFFF", size = 0.3) +
  scale_fill_gradientn(colors = rev(brewer.pal(11,'RdYlBu')), limits = c(-3.1,3.1)) +
  new_scale_fill() +
  geom_tile(data = scrna.df2, aes(x = scrna_cluster, y = 'cluster1', fill = cluster), height = 0.5, color = "#FFFFFF", size = 0.3) + 
  scale_fill_manual(values = clust_cols2$cols[order(as.character(clust_cols2$cluster))]) +
  new_scale_fill() +
  geom_tile(data = scrna.df2, aes(x = scrna_cluster, y = 'cluster2', fill = sort_id), height = 0.5, color = "#FFFFFF", size = 0.3) + 
  scale_fill_manual(values = sort_cols) +
  scale_x_discrete(limits = c('3_Il10neg','1_Il10neg','2_Il10neg','7_Il10neg','1_Il10pos','2_Il10pos','3_Il10pos','8_Il10neg','8_Il10pos','0_Il10pos','10_Il10pos','9_Il10neg','9_Il10pos','4_Il10neg','6_Il10neg','6_Il10pos','5_Il10neg','5_Il10pos')) +
  scale_y_discrete(limits = c('cluster1','cluster2','bI','bII','bIII','bIV','bV','bVI','bVII')) +
  plain.theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p2.scrna <- ggplot(data = scrna.df, aes(scrna_cluster, bulk_cluster)) +
  geom_tile(aes(fill = expr), color = "#FFFFFF", size = 0.3) +
  scale_fill_gradientn(colors = hcl.colors(256,'Viridis'), limits = c(.5,.9)) +
  scale_x_discrete(limits = c('3_Il10neg','1_Il10neg','2_Il10neg','1_Il10pos','2_Il10pos','3_Il10pos','8_Il10neg','8_Il10pos','0_Il10pos','10_Il10pos','9_Il10neg','9_Il10pos','4_Il10neg','6_Il10neg','6_Il10pos','5_Il10neg','5_Il10pos')) +
  scale_y_discrete(limits = c('Pseudo')) +
  plain.theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p.bulk <- ggplot(data = rna.df, aes(Celltype, bulk_cluster)) +
  geom_tile(aes(fill = expr), color = "#FFFFFF", size = 0.3) +
  scale_fill_gradientn(colors = rev(brewer.pal(11,'RdYlBu')), limits = c(-3.1,3.1)) +
  new_scale_fill() +
  geom_tile(data = rna.df2, aes(x = Celltype, y = 'cluster2', fill = sort_pop), height = 0.5, color = "#FFFFFF", size = 0.3) + 
  scale_fill_manual(values = bulk.cols[1:3]) +
  scale_x_discrete(limits = c('RNA_G1_neg','RNA_G2_neg','RNA_G3_neg','RNA_G1_pos','RNA_G2_pos','RNA_G3_pos','RNA_G1_stable','RNA_G2_stable','RNA_G3_stable')) +
  scale_y_discrete(limits = c('cluster1','cluster2','bI','bII','bIII','bIV','bV','bVI','bVII')) +
  plain.theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
plot_v1 <- p.bulk + p.scrna + plot_layout(ncol = 2, widths = c(2,4.3))
ggsave('comp_bulk_scrna_newgut.pdf', height = 6, width = 19, dpi = 300)

pa <- DimPlot(Il10gut, reduction = 'FDG', pt.size = 0.3, cols = clust_cols)
pb <- DimPlot(Il10gut, reduction = 'FDG', pt.size = 0.3, cols = sort_cols, group.by = 'SortID') + ylim(ggplot_build(pa)$layout$panel_scales_y[[1]]$range$range) + xlim(ggplot_build(pa)$layout$panel_scales_x[[1]]$range$range)
pc <- DimPlot(Il10gut, reduction = 'FDG', pt.size = 0.3, group.by = 'bulk_id', cols = bulk.cols)  + ylim(ggplot_build(pa)$layout$panel_scales_y[[1]]$range$range) + xlim(ggplot_build(pa)$layout$panel_scales_x[[1]]$range$range)
pd <- FeaturePlot(Il10gut,'entropy', reduction = 'FDG', pt.size = 0.3, order = TRUE) + scale_color_gradientn(colors = brewer.pal(9,'YlGnBu'))  + ylim(ggplot_build(pa)$layout$panel_scales_y[[1]]$range$range) + xlim(ggplot_build(pa)$layout$panel_scales_x[[1]]$range$range)
plot_v2 <- pa + pb + pc + pd
df1 <- FetchData(Il10gut, c('seurat_clusters','SortID','hash.ID','bulk_id','entropy','pseudo'))
pr.scrna <- ggplot(data = subset(df1, bulk_id %in% c('b_Il10neg','b_Il10pos','b_Il10stable')), aes(entropy,bulk_id)) +
  geom_density_ridges2(aes(fill = bulk_id), scale = 4) +
  scale_fill_manual(values = bulk.cols) +
  #scale_x_continuous(limits = c(0.0,0.1), expand = expansion(mult = c(0,0), add = c(0.01,0))) +
  plain.theme

Il10combined <- readRDS('Il10newcombined_ms.rds')
Il10combined$gut_cluster <- Il10gut$seurat_clusters[match(Cells(Il10combined), Cells(Il10gut), nomatch = NA)]
Il10combined$hash.IDv2 <- Il10combined$hash.ID
Il10combined$hash.IDv2[which(Cells(Il10combined) %in% Cells(Il10gut))] <- 'LILP.mesLN'
plot10 <- DimPlot(Il10combined, reduction = 'harmony', pt.size = 0.3, split.by = 'hash.IDv2', group.by = 'gut_cluster', cols = clust_cols) +
  theme(panel.border = element_rect(color = "#000000", fill = NA, size = 0.5), panel.grid = element_line(color = "#565656", linetype = 2, size = 0.2))
plot11 <- DimPlot(Il10combined, reduction = 'harmony', pt.size = 0.3, split.by = 'hash.IDv2', group.by = 'SortID', cols = sort_cols) +
  theme(panel.border = element_rect(color = "#000000", fill = NA, size = 0.5), panel.grid = element_line(color = "#565656", linetype = 2, size = 0.2))
plot12 <- plot10 + plot11 + plot_layout(nrow = 2)
ggsave(plot = plot12, filename = 'comp_lung_ms.pdf', height = 10, width = 13.5)
