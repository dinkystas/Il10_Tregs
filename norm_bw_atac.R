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
  as.numeric(read.delim(paste0('./ATAC_sam_out/' ,as.character(i),'.dupe_log.txt'), stringsAsFactors = FALSE, header = FALSE, sep = " ")[3,2])
}))
sample.info$scale_factor <- mean(sample.info$reads)/sample.info$reads
write.table('~/seq_programs/UCSC/faCount /media/dinky/Mendoza/mm39/mm39_reformed.fa','./ATAC_scripts/get_egs.sh', quote = FALSE, row.names = FALSE, col.names = FALSE)
system2('bash', args = './ATAC_scripts/get_egs.sh', stdout = '../../mm39/mm39_reformed_egs.txt')
egs <- read.delim('../../mm39/mm39_reformed_egs.txt', stringsAsFactors = FALSE)
egs <- sum(egs[nrow(egs),3:6])
final.commands <- unlist(lapply(1:nrow(sample.info),function(i){
  c1 <- paste0("bamCoverage -b ./ATAC_sam_out/", sample.info$Sample[i], ".duprm.bam --effectiveGenomeSize ", egs, " -bs 1 -p 4 --maxFragmentLength 131172 --scaleFactor ", sample.info$scale_factor[i], " -o ./ATAC_sam_out/", sample.info$Sample[i], ".normdt.bw")
}))
write.table(final.commands,'./ATAC_scripts/bam_bw.sh', quote = FALSE, row.names = FALSE, col.names = FALSE)

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
sample.info <- sample.info[-grep('RNA2', sample.info$Sample),]
sample.info$reads <- unlist(lapply(sample.info$Sample,function(i){
  as.numeric(read.delim(paste0('./RNA_sam_out/' ,as.character(i),'.dupe_log.txt'), stringsAsFactors = FALSE, header = FALSE, sep = " ")[3,2])
}))
sample.info$scale_factor <- mean(sample.info$reads)/sample.info$reads
final.commands <- unlist(lapply(1:nrow(sample.info),function(i){
  c1 <- paste0("bamCoverage -b ./RNA_sam_out/", sample.info$Sample[i], ".duprm.bam --effectiveGenomeSize ", egs, " -bs 1 -p 4 --scaleFactor ", sample.info$scale_factor[i], " -o ./RNA_sam_out/", sample.info$Sample[i], ".normdt.bw")
}))
write.table(final.commands,'./RNA_scripts/bam_bw.sh', quote = FALSE, row.names = FALSE, col.names = FALSE)