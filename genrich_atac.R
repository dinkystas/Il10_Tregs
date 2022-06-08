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
sample.info$reads <- unlist(lapply(sample.info$Sample,function(i){
  as.numeric(read.delim(paste0('./ATAC_sam_out/' ,as.character(i),'.dupe_log.txt'), stringsAsFactors = FALSE, header = FALSE, sep = " ")[3,2])/2
}))

auc.val <- round(tapply(sample.info$reads, sample.info$Celltype, mean)*0.00002)
dir.create('./Gen_out/', showWarnings = FALSE)
dir.create('./Gen_out/Peak', showWarnings = FALSE)
dir.create('./Gen_out/Log', showWarnings = FALSE)
cmd.1 <-  paste0('(Genrich -t ./ATAC_sam_out/', lapply(tapply(sample.info$Sample, sample.info$Celltype, paste0, '.byname.bam'), paste, collapse = ',./ATAC_sam_out/'),
                 ' -o ./Gen_out/Peak/', levels(sample.info$Celltype),
                 '.narrowPeak -j -d 25 -g 5 -v -q 0.01 -a ', auc.val, ') 2>./Gen_out/Log/', levels(sample.info$Celltype), '.stdout.txt')
cmd.2 <- paste('cat ./Gen_out/Peak/neg.narrowPeak ./Gen_out/Peak/pos.narrowPeak ./Gen_out/Peak/stable.narrowPeak > ./Gen_out/Peak/LI_comb.narrowPeak')
cmd.3 <- paste('sort -k1,1 -k2,2n ./Gen_out/Peak/LI_comb.narrowPeak > ./Gen_out/Peak/LI_comb_sort.narrowPeak')
cmd.4 <- paste('bedtools cluster -d -1 -i ./Gen_out/Peak/LI_comb_sort.narrowPeak > ./Gen_out/Peak/LI_clustered.narrowPeak')
write.table(c(cmd.1,cmd.2,cmd.3,cmd.4),paste0('./ATAC_scripts/Genrich.sh'), row.names = FALSE, col.names = FALSE, quote = FALSE)
