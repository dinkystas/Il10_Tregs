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
dir.create('./ATAC_scripts/', showWarnings = FALSE)
dir.create('./ATAC_TRIM_out/', showWarnings = FALSE)
dir.create('/media/dinky/Mendoza/mm39/mm39_100', showWarnings = FALSE)
cmd.gen <- "STAR --runThreadN 4 --runMode genomeGenerate --genomeDir /media/dinky/Mendoza/mm39/mm39_100 --genomeFastaFiles /media/dinky/Mendoza/mm39/mm39_reformed.fa"
cmd.0 <- paste0("TrimmomaticPE ",fastq.path, sample.info$Read1, " ",fastq.path, sample.info$Read2," -baseout ./ATAC_TRIM_out/", sample.info$Sample,".fastq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36")
dir.create('./ATAC_STAR_out', showWarnings = FALSE)
cmd.1 <- paste0('STAR --runThreadN 6 --runMode alignReads --genomeLoad NoSharedMemory --readFilesCommand zcat --genomeDir /media/dinky/Mendoza/mm39/mm39_100 --readFilesIn ./ATAC_TRIM_out/', sample.info$Sample,
                '_1P.fastq.gz ./ATAC_TRIM_out/', sample.info$Sample,
                '_2P.fastq.gz --outFileNamePrefix ./ATAC_STAR_out/', sample.info$Sample,
                ' --outSAMtype BAM Unsorted --outBAMcompression 6 --outFilterMultimapNmax 1 --outFilterMismatchNoverLmax 0.06 --outFilterMatchNminOverLread 0.35 --outFilterMatchNmin 30 --alignIntronMax 1 --alignEndsType Local')
cmd.1b <- paste0('rm ./ATAC_TRIM_out/', sample.info$Sample, '_1P.fastq.gz ./ATAC_TRIM_out/', sample.info$Sample, '_2P.fastq.gz ./ATAC_TRIM_out/', sample.info$Sample, '_1U.fastq.gz ./ATAC_TRIM_out/', sample.info$Sample, '_2U.fastq.gz')
dir.create('./ATAC_sam_out/', showWarnings = FALSE)
cmd.s0 <- paste0('samtools sort -@ 4 -n -o ./ATAC_sam_out/', sample.info$Sample, '.bam ./ATAC_STAR_out/', sample.info$Sample, 'Aligned.out.bam')
cmd.s1 <- paste0('rm ./ATAC_STAR_out/', sample.info$Sample, 'Aligned.out.bam')
cmd.s2 <- paste0('samtools fixmate -@ 4 -rm ./ATAC_sam_out/', sample.info$Sample, '.bam ./ATAC_sam_out/', sample.info$Sample, '.fixmate.bam')
cmd.s3 <- paste0('rm ./ATAC_sam_out/', sample.info$Sample, '.bam')
cmd.s4 <- paste0('samtools sort -@ 4 -o ./ATAC_sam_out/', sample.info$Sample, '.resort.bam ./ATAC_sam_out/', sample.info$Sample, '.fixmate.bam')
cmd.s5 <- paste0('rm ./ATAC_sam_out/', sample.info$Sample, '.fixmate.bam')
cmd.s6 <- paste0('(samtools markdup -@ 4 -l 1500 -r -d 100 -s ./ATAC_sam_out/', sample.info$Sample, '.resort.bam ./ATAC_sam_out/', sample.info$Sample, '.duprm.bam) 2>./ATAC_sam_out/', sample.info$Sample, '.dupe_log.txt')
cmd.s7 <- paste0('rm ./ATAC_sam_out/', sample.info$Sample, '.resort.bam')
cmd.s8 <- paste0('samtools sort -@ 4 -n -o ./ATAC_sam_out/', sample.info$Sample, '.byname.bam ./ATAC_sam_out/', sample.info$Sample, '.duprm.bam')
cmd.s9 <- paste0('samtools index -@ 4 -b ./ATAC_sam_out/', sample.info$Sample, '.duprm.bam')
write.table(c(cmd.gen,rbind(cmd.0,cmd.1,cmd.1b,cmd.s0,cmd.s1,cmd.s2,cmd.s3,cmd.s4,cmd.s5,cmd.s6,cmd.s7,cmd.s8,cmd.s9)),paste0('./ATAC_scripts/TRIM_STAR_samtools.sh'), row.names = FALSE, col.names = FALSE, quote = FALSE)

