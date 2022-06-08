library(stringr)
test.clustered <- data.frame(chr = read.delim('./Gen_out/Peak/LI_clustered.narrowPeak', header = FALSE, stringsAsFactors = FALSE)[,1],
                             start = read.delim('./Gen_out/Peak/LI_clustered.narrowPeak', header = FALSE, stringsAsFactors = FALSE)[,2],
                             summit = read.delim('./Gen_out/Peak/LI_clustered.narrowPeak', header = FALSE, stringsAsFactors = FALSE)[,2] + read.delim('./Gen_out/Peak/LI_clustered.narrowPeak', header = FALSE, stringsAsFactors = FALSE)[,10],
                             end = read.delim('./Gen_out/Peak/LI_clustered.narrowPeak', header = FALSE, stringsAsFactors = FALSE)[,3],
                             width = read.delim('./Gen_out/Peak/LI_clustered.narrowPeak', header = FALSE, stringsAsFactors = FALSE)[,3] - read.delim('./Gen_out/Peak/LI_clustered.narrowPeak', header = FALSE, stringsAsFactors = FALSE)[,2],
                             cluster = read.delim('./Gen_out/Peak/LI_clustered.narrowPeak', header = FALSE, stringsAsFactors = FALSE)[,11], stringsAsFactors = FALSE)
clusters <- data.frame(table(test.clustered$cluster))
colnames(clusters) <- c('cluster','size')
clusters$cluster <- as.integer(clusters$cluster)
testnumbs <- seq_along(clusters$size)
all.olaps <- lapply(testnumbs,function(h){
  curr <- test.clustered[which(test.clustered$cluster == h),]
  curr$peak <- paste0('cluster_', h, '_peak_',seq_along(curr$width))
  all.olaps <- unlist(lapply(seq_along(curr$peak),function(i){
    summit.a <- curr$summit[i]
    length(which(curr$start < summit.a & curr$end > summit.a))/length(curr$peak)
  }))
})
clusters$olap.frac <- unlist(lapply(all.olaps,function(i){
  length(which(i == 1))/length(i)
}))
atlas.final <- do.call(rbind,lapply(testnumbs,function(h){
  curr <- test.clustered[which(test.clustered$cluster == h),]
  if(clusters$olap.frac[h] == 1){
    data.frame(chr = curr$chr[1], start = round(mean(curr$start)), summit = round(mean(curr$summit)), end = round(mean(curr$end)), width = (round(mean(curr$end)) - round(mean(curr$start))))
  }
  else{
    curr <- test.clustered[which(test.clustered$cluster == h),]
    all.coords <- data.frame(rbind(cbind(curr$start,'start'),cbind(curr$summit,'summit'),cbind(curr$end,'end')), stringsAsFactors = FALSE)
    all.coords$X1 <- as.numeric(all.coords$X1)
    all.coords <- all.coords[order(all.coords[,1]),]
    summits <- which(all.coords$X2 == 'summit')
    test <- do.call(rbind,lapply(summits,function(i){
      before.pos <- (i-1)
      after.pos <- (i+1)
      summits.tomean <- as.integer(i)
      while(all.coords[before.pos,2] == 'summit'){
        summits.tomean <- c(summits.tomean,before.pos)
        before.pos <- (before.pos - 1)}
      start.temp <- all.coords[before.pos,1]
      while(all.coords[after.pos,2] == 'summit'){
        summits.tomean <- c(summits.tomean,after.pos)
        after.pos <- (after.pos + 1)}
      end.temp <- all.coords[after.pos,1]
      summit.temp <- round(mean(all.coords[summits.tomean,1]))
      data.frame(chr = curr$chr[1], start = start.temp, summit = summit.temp, end = end.temp, stringsAsFactors = FALSE)
    }))
    test <- unique(test)
    test$width <- test$end - test$start
    test
  }
}))
atlas.final <- atlas.final[-grep("_", atlas.final$chr),]
atlas.final <- atlas.final[-which(atlas.final$width > 3500 & atlas.final$width < 16000),]
atlas.final$chr[which(atlas.final$chr == 'chrX')] <- 'chr20'
atlas.final$chr[which(atlas.final$chr == 'chrY')] <- 'chr21'
atlas.final$chr[which(atlas.final$chr == 'chrM')] <- 'chr22'
atlas.final$chr <- unlist(lapply(1:nrow(atlas.final),function(i){as.numeric(substring(atlas.final$chr[i],4,5))}))
atlas.final <- atlas.final[order(atlas.final$chr, atlas.final$start),]
atlas.final$chr <- paste0('chr',atlas.final$chr)
atlas.final$chr[which(atlas.final$chr == 'chr20')] <- 'chrX'
atlas.final$chr[which(atlas.final$chr == 'chr21')] <- 'chrY'
atlas.final$chr[which(atlas.final$chr == 'chr22')] <- 'chrM'
rownames(atlas.final) <- 1:nrow(atlas.final)
atlas.final$peak <- paste0('peak_',seq_along(atlas.final$width))
atlas.final$chr <- as.character(atlas.final$chr)
atlas.all <- data.frame(chr = atlas.final$chr, start = atlas.final$start, end = atlas.final$end, peak = atlas.final$peak)
write.table(atlas.all, './Gen_out/Peak/R_LI_atlas_all.bed', sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table('bedtools getfasta -fo ./Gen_out/Peak/R_LI_atlas_all_seqs.txt -tab -fi /media/dinky/Mendoza/mm39/mm39_reformed.fa -bed ./Gen_out/Peak/R_LI_atlas_all.bed','./ATAC_scripts/get_peak_seqs.sh', quote = FALSE, row.names = FALSE, col.names = FALSE)
system2('bash', args = './ATAC_scripts/get_peak_seqs.sh')
peak.seqs <- read.delim('./Gen_out/Peak/R_LI_atlas_all_seqs.txt', stringsAsFactors = FALSE, header = FALSE)
peak.seqs$replen <- unlist(lapply(1:nrow(peak.seqs),function(i){
  str_count(peak.seqs[i,2],"[a-z]")/nchar(peak.seqs[i,2])
}))
tmp.peaks <- data.frame(peak = atlas.all$peak, loc = paste0(atlas.all$chr,":",atlas.all$start,"-",atlas.all$end), stringsAsFactors = FALSE)
peak.seqs <- merge(tmp.peaks, peak.seqs, by.x = 'loc', by.y = 'V1')
peak.seqs <- subset(peak.seqs, replen < 0.7)
write.table(peak.seqs, './Gen_out/Peak/R_LI_atlas_norep_seqs.txt', sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
rownames(atlas.final) <- atlas.final$peak
atlas.final <- atlas.final[peak.seqs$peak,]
rm(all.olaps,atlas.all,clusters,peak.seqs,test.clustered,testnumbs,tmp.peaks)
library(GenomicAlignments)
library(GenomicRanges)
library(GenomicFeatures)
library(ChIPpeakAnno)
mouse.genes.anno <- readRDS('mm39_gtf.rds')
mouse.genes.anno <- genes(makeTxDbFromGRanges(mouse.genes.anno))
atlas.GR <- atlas.final[,c(1,2,4,6)]
colnames(atlas.GR) <- c('chrom','start','end','name')
atlas.GR <- toGRanges(data = atlas.GR, format = 'BED')
annot.aprom <- annoPeaks(atlas.GR, annoData=mouse.genes.anno, select = 'bestOne', bindingType = "startSite", bindingRegion = c(-2000,500))
annot.aprom <- data.frame(peak = annot.aprom@elementMetadata$peak, Gene = annot.aprom@elementMetadata$feature, type = 'promoter', stringsAsFactors = FALSE)
annot.aintra <- annoPeaks(atlas.GR, annoData=mouse.genes.anno, select = 'bestOne', bindingType = "fullRange", bindingRegion = c(0,1))
annot.aintra <- data.frame(peak = annot.aintra@elementMetadata$peak, Gene = annot.aintra@elementMetadata$feature, type = 'intragenic', stringsAsFactors = FALSE)
annot.ainter <- annoPeaks(atlas.GR, annoData=mouse.genes.anno, select = 'bestOne', bindingType = "fullRange", bindingRegion = c(-100000,100000))
annot.ainter <- data.frame(peak = annot.ainter@elementMetadata$peak, Gene = annot.ainter@elementMetadata$feature, type = 'intergenic', stringsAsFactors = FALSE)
testype <- data.frame(peak_og = atlas.final$peak, annot.aprom[match(atlas.final$peak, annot.aprom$peak),], stringsAsFactors = FALSE)
testype[which(is.na(testype$peak)),] <- cbind(testype$peak_og[which(is.na(testype$peak))],annot.aintra[match(testype$peak_og[which(is.na(testype$peak))], annot.aintra$peak),])
testype[which(is.na(testype$peak)),] <- cbind(testype$peak_og[which(is.na(testype$peak))],annot.ainter[match(testype$peak_og[which(is.na(testype$peak))], annot.ainter$peak),])
testype[which(is.na(testype$peak)),] <- cbind(testype$peak_og[which(is.na(testype$peak))],testype$peak_og[which(is.na(testype$peak))], rep(NA,length(testype$peak_og[which(is.na(testype$peak))])), rep('None',length(testype$peak_og[which(is.na(testype$peak))])))
atlas.final <- merge(atlas.final, testype[,2:4], by = "peak", sort = FALSE)
atlas.final$chr[which(atlas.final$chr == 'chrX')] <- 'chr20'
atlas.final$chr[which(atlas.final$chr == 'chrY')] <- 'chr21'
atlas.final$chr[which(atlas.final$chr == 'chrM')] <- 'chr22'
atlas.final$chr <- unlist(lapply(1:nrow(atlas.final),function(i){
  as.numeric(substring(atlas.final$chr[i],4,5))
}))
atlas.final <- atlas.final[order(atlas.final$chr, atlas.final$start),]
atlas.final$chr <- paste0('chr',atlas.final$chr)
atlas.final$chr[which(atlas.final$chr == 'chr20')] <- 'chrX'
atlas.final$chr[which(atlas.final$chr == 'chr21')] <- 'chrY'
atlas.final$chr[which(atlas.final$chr == 'chr22')] <- 'chrM'
write.table(atlas.final[,c(2,3,5,1)], './Gen_out/Peak/R_LI_atlas_formal.bed', sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(atlas.final, './Gen_out/Peak/R_LI_atlas_winfo.bed', sep = "\t", row.names = FALSE, quote = FALSE)
