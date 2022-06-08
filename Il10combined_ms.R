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
Il10combined <- ScaleData(Il10combined, features = rownames(Il10combined))
Il10combined <- FindVariableFeatures(Il10combined, selection.method = "vst", nfeatures = 3000)
Il10combined <- RunPCA(Il10combined, features = VariableFeatures(object = Il10combined))
Il10combined <- FindNeighbors(Il10combined, reduction = 'pca', dims = 1:30)
Il10combined <- FindClusters(Il10combined, resolution = 0.5)
Il10combined <- RunUMAP(Il10combined, reduction = 'pca', dims = 1:30)
df1 <- runFDG(Il10combined@reductions$pca@cell.embeddings[,1:30],Il10combined@graphs$RNA_snn, python.addr = "python")
colnames(df1) <- c('FDG_1','FDG_2')
Il10combined[['FDG']] <- CreateDimReducObject(as.matrix(df1))
rm(df1)
og.markers <- FindAllMarkers(Il10combined, logfc.threshold = 1)
og.markers <- og.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
t1 <- as.matrix(table(Il10combined$seurat_clusters,Il10combined$hash.ID))
t1 <- intersect(rownames(as.data.frame.matrix(t1/rowSums(t1)))[which(as.data.frame.matrix(t1/rowSums(t1))$mesLN > .8)],
                unique(as.character(og.markers$cluster[grep('Bcl2',og.markers$gene)])))
t1 <- Cells(subset(Il10combined, seurat_clusters == t1))
Il10combined <- RunHP(object = Il10combined, early.cell = sample(t1,1), n.iter = 600, nPCs = 30)
saveRDS(Il10combined,'Il10newcombined_ms.rds')
