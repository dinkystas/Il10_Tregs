library(GenomicAlignments)
library(GenomicRanges)
library(GenomicFeatures)
library(ChIPpeakAnno)
library(DESeq2)
fastq.path <- "../Seq_files/ATAC_FASTQ/"
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
sample.info$reads <- unlist(lapply(sample.info$Sample,function(i){
  as.numeric(read.delim(paste0('./ATAC_sam_out/' ,as.character(i),'.dupe_log.txt'), stringsAsFactors = FALSE, header = FALSE, sep = " ")[3,2])/2
}))
sample.info$Bamfile <- paste0('./ATAC_sam_out/',sample.info$Sample,'.duprm.bam')
atlas.GR <- toGRanges(data = './Gen_out/Peak/R_LI_atlas_formal.bed', format = 'BED')
atlas.final <- read.delim('./Gen_out/Peak/R_LI_atlas_winfo.bed', stringsAsFactors = FALSE)
#peak.counts <- summarizeOverlaps(atlas.GR, sample.info$Bamfile)
#saveRDS(peak.counts, 'atac_sum_olap.rds')
peak.counts <- readRDS('atac_sum_olap.rds')
rawcounts <- as.data.frame(peak.counts@assays@data@listData$counts)
rownames(rawcounts) <- peak.counts@rowRanges@ranges@NAMES
sample.info$reads_chrM <- unlist(rawcounts[atlas.final$peak[which(atlas.final$chr == 'chrM')],])
peak.counts <- peak.counts[-which(peak.counts@rowRanges@ranges@NAMES == 'peak_77352'),]
atlas.final <- atlas.final[-which(atlas.final$peak == 'peak_77352'),]
sample.info$reads_align <- colSums(rawcounts)
saveRDS(rawcounts, 'atac_raw_counts.rds')
colnames(peak.counts) <- sample.info$Sample
colData(peak.counts)$Celltype <- sample.info$Celltype
colData(peak.counts)$Replicate <- sample.info$Replicate
peak.counts <- DESeqDataSet(peak.counts, design = ~ Celltype + Replicate)
#sizeFactors(peak.counts) <- sample.info$reads_chrM/mean(sample.info$reads_chrM)
fpkm.counts <- fpkm(peak.counts, robust = TRUE)
saveRDS(fpkm.counts,paste0('atac_fpkm.rds'))
fpkm.counts <- fpkm(peak.counts, robust = FALSE)
saveRDS(fpkm.counts,paste0('atac_fpkm_raw.rds'))
peak.counts <- DESeq(peak.counts)
sample.info$size_factors <- sizeFactors(peak.counts)
#write.table('~/seq_programs/UCSC/faCount /media/dinky/Mendoza/mm39/mm39_reformed.fa','./ATAC_scripts/get_egs.sh', quote = FALSE, row.names = FALSE, col.names = FALSE)
#system2('bash', args = './ATAC_scripts/get_egs.sh', stdout = '../../mm39/mm39_reformed_egs.txt')
egs <- read.delim('../../mm39/mm39_reformed_egs.txt', stringsAsFactors = FALSE)
egs <- sum(egs[nrow(egs),3:6])
final.commands <- unlist(lapply(1:nrow(sample.info),function(i){
  c1 <- paste0("bamCoverage -b ./ATAC_sam_out/", sample.info$Sample[i], ".duprm.bam --effectiveGenomeSize ", egs, " -bs 1 -p 4 --maxFragmentLength 131172 --scaleFactor ", 1/sample.info$size_factors[i], " -o ./ATAC_sam_out/", sample.info$Sample[i], ".normdt.bw")
}))
write.table(final.commands,'./ATAC_scripts/bam_bw.sh', quote = FALSE, row.names = FALSE, col.names = FALSE)
comp.combn <- data.frame(combn(rev(unique(as.character(peak.counts$Celltype))),2))
atac.results <- lapply(comp.combn, function(i){
  temp.res <- as.data.frame(results(peak.counts, c('Celltype',as.character(i[1]),as.character(i[2]))))
  temp.res <- data.frame(peak = rownames(temp.res), temp.res)
  temp.res <- merge(atlas.final, temp.res, by.x = "peak", by.y = "peak", sort = FALSE)
  temp.res <- temp.res[,c(1:ncol(atlas.final),(ncol(atlas.final) + 2),ncol(temp.res))]
  colnames(temp.res) <- c('peak','chr','start','summit','end','width','Gene','type',paste0('lg2FC_',substring(i[1],1,1),'v',substring(i[2],1,1)),paste0('padj_',substring(i[1],1,1),'v',substring(i[2],1,1)))
  return(temp.res)
})
atac.results <- do.call('cbind',atac.results)
atac.results <- atac.results[,!duplicated(substring(colnames(atac.results),regexpr('\\.', colnames(atac.results))+1,nchar(colnames(atac.results))))]
colnames(atac.results) <- substring(colnames(atac.results),regexpr('\\.', colnames(atac.results))+1,nchar(colnames(atac.results)))
normcounts <- counts(peak.counts, normalized = TRUE)
normcounts <- log2(normcounts + 1)
testaves <- do.call(rbind,lapply(seq_along(rownames(normcounts)),function(i){
  tapply(as.vector(normcounts[i,]),sample.info$Celltype,mean)
}))
colnames(testaves) <- c('neg_lg2tag','pos_lg2tag','stable_lg2tag')
rownames(testaves) <- rownames(normcounts)
atac.results <- data.frame(atac.results[,1:8],testaves,atac.results[,9:14])
saveRDS(atac.results,paste0('atac_results.rds'))
saveRDS(normcounts,paste0('atac_normcounts.rds'))
