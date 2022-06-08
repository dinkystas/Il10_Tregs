library(GenomicAlignments)
library(GenomicRanges)
library(GenomicFeatures)
library(ChIPpeakAnno)
library(DESeq2)
library(ggplot2)
library(patchwork)
plain.theme <- theme(axis.line.x = element_line(color = '#000000', size = 1, lineend = 'round'),
                     axis.line.y = element_line(color = '#000000', size = 1, lineend = 'round'),
                     axis.ticks = element_line(color = 'black', size = 1, lineend = 'round'),
                     panel.background = element_rect(fill = 'white'),
                     text = element_text(size = 10, face = 'bold', color = 'black'),
                     axis.text = element_text(size = 10, face = 'bold', color = 'black'))
fastq.path <- "../Seq_files/RNA_FASTQ/"
all.fastqs <- list.files(path = fastq.path, pattern = '.fastq.gz')
Lib.name <- unique(unlist(lapply(seq_along(all.fastqs),function(i){
  paste0(substring(all.fastqs[i],1,gregexpr('IGO', all.fastqs[i])[[1]][1]-2))
})))
Lib.R1 <- all.fastqs[grep('R1', all.fastqs)]
Lib.R2 <- all.fastqs[grep('R2', all.fastqs)]
sample.info <- data.frame(Sample = Lib.name, Read1 = Lib.R1, Read2 = Lib.R2, stringsAsFactors = FALSE)
sample.info$Replicate <- unlist(lapply(seq_along(sample.info$Sample),function(i){
  substring(sample.info$Sample[i],gregexpr('_',sample.info$Sample[i])[[1]][1]+1,gregexpr('_',sample.info$Sample[i])[[1]][2]-1)
}))
sample.info$Celltype <- unlist(lapply(seq_along(sample.info$Sample),function(i){
  substring(sample.info$Sample[i],gregexpr('_',sample.info$Sample[i])[[1]][2]+1,1000)
}))
sample.info$Celltype <- factor(sample.info$Celltype)
sample.info$Replicate <- factor(sample.info$Replicate)
sample.info <- sample.info[-grep('RNA2', sample.info$Sample),]
sample.info$reads <- unlist(lapply(sample.info$Sample,function(i){
  as.numeric(read.delim(paste0('./RNA_sam_out/' ,as.character(i),'.dupe_log.txt'), stringsAsFactors = FALSE, header = FALSE, sep = " ")[3,2])
}))
sample.info$Bamfile <- paste0('./RNA_sam_out/',sample.info$Sample,'.duprm.bam')
mouse.genes.anno <- readRDS('mm39_gtf.rds')
mouse.genes.anno <- genes(makeTxDbFromGRanges(mouse.genes.anno))
#rna.counts <- summarizeOverlaps(mouse.genes.anno, sample.info$Bamfile)
#saveRDS(rna.counts, 'rna_sum_olap.rds')
rna.counts <- readRDS('rna_sum_olap.rds')
rawcounts <- as.data.frame(rna.counts@assays@data@listData$counts)
rownames(rawcounts) <- rna.counts@rowRanges@ranges@NAMES
sample.info$reads_mt <- colSums(rawcounts[grep('mt-Rnr',rownames(rawcounts)),])
sample.info$reads_align <- colSums(rawcounts)
#saveRDS(rawcounts, 'rna_raw_counts.rds')
colnames(rna.counts) <- sample.info$Sample
colData(rna.counts)$Celltype <- sample.info$Celltype
colData(rna.counts)$Replicate <- sample.info$Replicate
rna.counts <- DESeqDataSet(rna.counts, design = ~ Celltype + Replicate)
#sizeFactors(rna.counts) <- sample.info$reads_mt/mean(sample.info$reads_mt)
fpkm.counts <- fpkm(rna.counts, robust = TRUE)
rna.counts <- DESeq(rna.counts)
sample.info$size_factors <- sizeFactors(rna.counts)
egs <- read.delim('../../mm39/mm39_reformed_egs.txt', stringsAsFactors = FALSE)
egs <- sum(egs[nrow(egs),3:6])
final.commands <- unlist(lapply(1:nrow(sample.info),function(i){
  c1 <- paste0("bamCoverage -b ./RNA_sam_out/", sample.info$Sample[i], ".duprm.bam --effectiveGenomeSize ", egs, " -bs 1 -p 4 --scaleFactor ", 1/sample.info$size_factors[i], " -o ./RNA_sam_out/", sample.info$Sample[i], ".normdt.bw")
}))
write.table(final.commands,'./RNA_scripts/bam_bw.sh', quote = FALSE, row.names = FALSE, col.names = FALSE)
comp.combn <- data.frame(combn(rev(unique(as.character(rna.counts$Celltype))),2))
rna.results <- lapply(comp.combn, function(i){
  temp.res <- as.data.frame(results(rna.counts, c('Celltype',as.character(i[1]),as.character(i[2]))))
  temp.res <- data.frame(Gene = rownames(temp.res), temp.res[,c(2,6)])
  colnames(temp.res) <- c('Gene',paste0('lg2FC_',substring(i[1],1,1),'v',substring(i[2],1,1)),paste0('padj_',substring(i[1],1,1),'v',substring(i[2],1,1)))
  return(temp.res)
})
rna.results <- do.call('cbind',rna.results)
rna.results <- rna.results[,!duplicated(substring(colnames(rna.results),regexpr('\\.', colnames(rna.results))+1,nchar(colnames(rna.results))))]
colnames(rna.results) <- substring(colnames(rna.results),regexpr('\\.', colnames(rna.results))+1,nchar(colnames(rna.results)))
normcounts <- counts(rna.counts, normalized = TRUE)
normcounts <- log2(normcounts + 1)
testaves <- do.call(rbind,lapply(seq_along(rownames(normcounts)),function(i){
  tapply(as.vector(normcounts[i,]),sample.info$Celltype,mean)
}))
colnames(testaves) <- c('neg_lg2expr','pos_lg2expr','stable_lg2expr')
rownames(testaves) <- rownames(normcounts)
rna.results <- data.frame(Gene = rna.results[,1],testaves,rna.results[,2:7], stringsAsFactors = FALSE)
saveRDS(fpkm.counts,paste0('rna_fpkm.rds'))
saveRDS(normcounts,paste0('rna_normcounts.rds'))
saveRDS(rna.results,paste0('rna_results.rds'))