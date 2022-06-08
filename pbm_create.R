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
rm(all.fastqs,fastq.path,Lib.name,Lib.R1,Lib.R2)
elbplot <- function(vals){
  tempdf <- data.frame(torank = as.vector(vals))
  tempdf <- data.frame(torank = tempdf[order(tempdf$torank, decreasing = TRUE),])
  tempdf$Rank <- c(1:nrow(tempdf))
  tempdf <- subset(tempdf, torank > 0)
  plot(tempdf$Rank,tempdf$torank)
}
fpkm.counts <- readRDS(paste0('rna_fpkm.rds'))
rna.results <- readRDS(paste0('rna_results.rds'))
grp.fpkm <- do.call(rbind,lapply(1:nrow(fpkm.counts),function(i){
  tapply(fpkm.counts[i,], sample.info$Celltype,mean, na.rm = TRUE)
}))
rownames(grp.fpkm) <- rownames(fpkm.counts)
toexclude <- c('Ighv','Ighj','Ighd','Igkv','Igkj','Iglv','Iglj','Trav','Traj','Trbv','Trbj','Trbd','Trgv','Trgj','Trdv','Trdj','Trdd')
toexclude <- unique(c(unlist(lapply(toexclude,function(i){
  rownames(rna.results)[grep(i,rownames(rna.results))]
}))))
toexclude <- unique(c(toexclude,rownames(fpkm.counts)[which(rowMeans(fpkm.counts) <= mean(fpkm.counts['Cd8a',]))]))
rna.results <- rna.results[-match(toexclude,rownames(rna.results)),]
fpkm.counts <- fpkm.counts[-match(toexclude,rownames(fpkm.counts)),]
grp.fpkm <- grp.fpkm[-match(toexclude,rownames(grp.fpkm)),]
atac.results <- readRDS('atac_results.rds')
grp.co <- unique(c(rownames(grp.fpkm)[which(grp.fpkm[,1] > 1)],rownames(grp.fpkm)[which(grp.fpkm[,2] > 1)],rownames(grp.fpkm)[which(grp.fpkm[,3] > 1)]))
tf.info <- read.delim('../../mm39/TF_Information.txt', stringsAsFactors = FALSE)[,c(4,7,8,9,10)]
tf.info <- tf.info[which(tf.info$TF_Name %in% rna.results$Gene),]
tf.info <- tf.info[grep('M', tf.info$Motif_ID),]
grp.co <- intersect(unique(tf.info$TF_Name),grp.co)
tf.info <- tf.info[which(tf.info$TF_Name %in% grp.co),]
tf.info$Unique_name <- make.unique(tf.info$TF_Name)
for.motifs <- do.call(rbind,lapply(seq_along(tf.info$Motif_ID),function(i){
  currnumb <- i
  temp.motif <- read.delim(paste0('../../mm39/pwms_all_motifs/', tf.info$Motif_ID[currnumb], '.txt'), row.names = 1)
  if(nrow(temp.motif) > 0){
    tfname <- tf.info$Unique_name[i]
    temp.motif <- rbind(c(paste0('>',paste0(tfname)),'','',''),unname(temp.motif[,1:4]))
    colnames(temp.motif) <- c('A','B','C','D')
    temp.motif}
  else{c(NA,NA,NA,NA)}
}))
for.motifs <- for.motifs[-which(is.na(for.motifs$A)),]
#dir.create('./meme_stuff', showWarnings = FALSE)
#write.table(for.motifs, file = "./meme_stuff/expr_tfs_final.motifs", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
peak.seqs <- read.delim('./Gen_out/Peak/R_LI_atlas_norep_seqs.txt', stringsAsFactors = FALSE, header = FALSE)
peak.seqs <- peak.seqs[-grep('peak_77352', peak.seqs$V2),]
peak.seqs.fa <- do.call(rbind,lapply(1:nrow(peak.seqs),function(i){
  rbind(paste0('>',peak.seqs$V2[i]),peak.seqs$V3[i])
}))
#write.table(peak.seqs.fa, './Gen_out/Peak/peak_seqs.fa', sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
rm(peak.seqs,peak.seqs.fa)
# more ./meme_stuff/expr_tfs_final.motifs | chen2meme > ./meme_stuff/expr_tfs_final.meme
# ame --verbose 2 --oc ./meme_stuff/ame --evalue-report-threshold 319 ./Gen_out/Peak/peak_seqs.fa ./meme_stuff/expr_tfs_final.meme
alltfs <- read.delim('./meme_stuff/ame/ame.tsv', stringsAsFactors = FALSE, header = TRUE)
alltfs.peaks <- as.numeric(na.omit(alltfs$pos))
alltfs.bkgd <- as.numeric(na.omit(alltfs$neg))
alltfs <- alltfs[1:(nrow(alltfs)-3),c(3,5,7,8,9,13:17)]
testa <- read.delim('./meme_stuff/ame/ame.tsv', header = TRUE, as.is = TRUE, colClasses = "character")[1:nrow(alltfs),7]
colnames(alltfs) <- c('motif_id','consensus','padj','eval','tests','pwm_min','peak_hits','tp_percent','bkg_hits','fp_percent')
alltfs$padj <- unlist(lapply(1:length(testa),function(i){
  if(grepl('e',testa[i])){as.numeric(substring(testa[i],regexpr('e',testa[i])+1,10000)) + log10(as.numeric(substring(testa[i],1,regexpr('e',testa[i])-1)))}
  else{
    if(grepl('E',testa[i])){as.numeric(substring(testa[i],regexpr('E',testa[i])+1,10000)) + log10(as.numeric(substring(testa[i],1,regexpr('E',testa[i])-1)))}
    else{log10(as.numeric(testa[i]))}}
}))
rm(testa)
alltfs$tp_percent <- alltfs$peak_hits/alltfs.peaks
alltfs$fp_percent <- alltfs$peak_hits/alltfs.bkgd
alltfs$tf_gene <- unlist(lapply(seq_along(alltfs$motif_id),function(i){
  subpos <- regexpr("\\.",alltfs$motif_id[i])[1]
  if(subpos == -1){substring(alltfs$motif_id[i],1,1000)}
  else{
    substring(alltfs$motif_id[i],1,subpos-1)}
}))
alltfs$motif_len <- nchar(alltfs$consensus)
rm(alltfs.bkgd,alltfs.peaks,toexclude)
tf.info <- tf.info[which(tf.info$Unique_name %in% alltfs$motif_id),]
temp1 <- tapply(alltfs$tp_percent,alltfs$tf_gene,max)
alltfs$best_percent <- unlist(lapply(1:nrow(alltfs),function(i){alltfs$tp_percent[i] == temp1[alltfs$tf_gene[i]]}))
alltfs <- alltfs[which(alltfs$best_percent),]

for.motifs2 <- do.call(rbind,lapply(seq_along(alltfs$motif_id),function(i){
  currnumb <- i
  currid <- tf.info$Motif_ID[which(tf.info$Unique_name == alltfs$motif_id[currnumb])]
  temp.motif <- read.delim(paste0('../../mm39/pwms_all_motifs/', currid, '.txt'), row.names = 1)
  tfname <- alltfs$tf_gene[i]
  temp.motif <- rbind(c(paste0('>',tfname),'','',''),unname(temp.motif[,1:4]))
  colnames(temp.motif) <- c('A','B','C','D')
  temp.motif
}))
#write.table(for.motifs2, file = "./meme_stuff/for_tomtom.motifs", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
# more ./meme_stuff/for_tomtom.motifs | chen2meme > ./meme_stuff/for_tomtom.meme
# tomtom -oc ./meme_stuff/tomtom -thresh 1 ./meme_stuff/for_tomtom.meme ./meme_stuff/for_tomtom.meme
tomtom.res <- read.delim('./meme_stuff/tomtom/tomtom.tsv', sep = "\t", stringsAsFactors = FALSE)
tomtom.res <- tomtom.res[1:(nrow(tomtom.res)-3),]
tomtom.res[grep('Nkx', tomtom.res$Query_ID),'Query_ID'] <- 'Nkx6-2'
tomtom.res[grep('Nkx', tomtom.res$Target_ID),'Target_ID'] <- 'Nkx6-2'
to.rm <- na.omit(unique(unlist(lapply(1:nrow(tomtom.res),function(i){
  if(tomtom.res$Target_ID[i] == tomtom.res$Query_ID[i]){return(i)}
  else{return(NA)}
}))))
tomtom.res$p_olap <- tomtom.res$Overlap/nchar(tomtom.res$Query_consensus)
tomtom.res <- tomtom.res[-to.rm,]
tt_fam <- lapply(seq_along(alltfs$tf_gene),function(i){
  cutoff <- .00001
  currtf <- alltfs$tf_gene[i]
  currfam <- unique(tf.info$Family_Name[which(tf.info$TF_Name == currtf)])
  currfam <- unique(tf.info$TF_Name[which(tf.info$Family_Name == currfam)])
  hits <- tomtom.res$Target_ID[which(tomtom.res$Query_ID == currtf & tomtom.res$E.value < cutoff)]
  finalhits <- unique(sort(intersect(currfam,c(currtf,hits))))
})
names(tt_fam) <- alltfs$tf_gene
prevlen <- length(unique(tt_fam))
currlen <- 0
while(currlen != prevlen){
  prevlen <- length(unique(tt_fam))
  tt_fam <- lapply(seq_along(tt_fam),function(h){
    currfam <- tt_fam[[h]]
    sort(unique(unlist(lapply(currfam,function(i){
      currtf <- as.character(i)
      unique(sort(unlist(tt_fam[which(unlist(lapply(tt_fam,function(j){currtf %in% j})))])))
    }))))
  })
  names(tt_fam) <- alltfs$tf_gene
  currlen <- length(unique(tt_fam))
}
alltfs$tt_fam <- unlist(lapply(1:nrow(alltfs),function(i){
  paste(tt_fam[[i]], collapse = ",")
}))
famkey <- unique(data.frame(fam_members = alltfs$tt_fam))
famkey$fam_id <- paste('family',1:nrow(famkey), sep = "_")
alltfs$fam_id <- unlist(lapply(1:nrow(alltfs),function(i){
  famkey$fam_id[which(famkey$fam_members == alltfs$tt_fam[i])]
}))
temp2 <- tapply(alltfs$tp_percent,alltfs$consensus,max)
alltfs$best_percent_2 <- unlist(lapply(1:nrow(alltfs),function(i){alltfs$tp_percent[i] == temp2[alltfs$consensus[i]]}))
alltfs.all <- alltfs
alltfs <- alltfs[which(alltfs$best_percent_2),]
alltfs <- alltfs[order(alltfs$tp_percent, decreasing = TRUE),]
alltfs <- alltfs[-which(duplicated(alltfs$consensus)),]
alltfs <- subset(alltfs, padj < log10(0.01))
rownames(alltfs) <- alltfs$tf_gene
currlen <- length(unique(alltfs$tt_fam))
rm(tomtom.res,temp1,temp2,to.rm)
for.motifs3 <- do.call(rbind,lapply(1:nrow(alltfs),function(i){
  currnumb <- i
  currid <- tf.info$Motif_ID[which(tf.info$Unique_name == alltfs$motif_id[currnumb])]
  temp.motif <- read.delim(paste0('../../mm39/pwms_all_motifs/', currid, '.txt'), row.names = 1)
  tfname <- alltfs$tf_gene[i]
  temp.motif <- rbind(c(paste0('>',tfname),'','',''),unname(temp.motif[,1:4]))
  colnames(temp.motif) <- c('A','B','C','D')
  temp.motif
}))
#write.table(for.motifs3, file = "./meme_stuff/post_tomtom.motifs", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
# more ./meme_stuff/post_tomtom.motifs | chen2meme > ./meme_stuff/post_tomtom.meme
# fimo --max-stored-scores 10000000  --oc ./meme_stuff/fimo --thresh 0.0001 --verbosity 2 --max-strand ./meme_stuff/post_tomtom.meme ./Gen_out/Peak/peak_seqs.fa
first.fimo.res <- read.delim('./meme_stuff/fimo/fimo.tsv', sep = "\t", stringsAsFactors = FALSE)
first.fimo.res <- first.fimo.res[1:(nrow(first.fimo.res)-3),c(1,3,4,5,7,8,9,10)]
first.fimo.res <- first.fimo.res[which(first.fimo.res$sequence_name %in% atac.results$peak),]
first.fimo.res$motif_len <- nchar(first.fimo.res$matched_sequence)
first.fimo.res$norm_score <- 1
alltfs <- alltfs[which(alltfs$tf_gene %in% unique(first.fimo.res$motif_id)),]
currlen <- length(unique(alltfs$tt_fam))
first.fimo.res$family <- alltfs$fam_id[match(first.fimo.res$motif_id, alltfs$tf_gene)]
test2 <- as.matrix(data.frame(cbind(table(first.fimo.res$family,first.fimo.res$sequence_name))))
alltfs$percent <- unlist(lapply(1:nrow(alltfs),function(i){
  fam <- alltfs$fam_id[i]
  length(which(test2[fam,] > 0L))/ncol(test2)
}))
test3 <- as.matrix(data.frame(cbind(table(first.fimo.res$motif_id,first.fimo.res$sequence_name))))
alltfs$percent_motif <- unlist(lapply(1:nrow(alltfs),function(i){
  fam <- alltfs$tf_gene[i]
  length(which(test3[fam,] > 0L))/ncol(test3)
}))
alltfs <- alltfs[which(alltfs$percent > 0.02),]
currlen <- length(unique(alltfs$tt_fam))
first.fimo.res <- first.fimo.res[which(first.fimo.res$motif_id %in% alltfs$motif_id),]
first.fimo.res$pbm_id <- interaction(first.fimo.res$family,first.fimo.res$sequence_name)
pbm <- as.vector(tapply(first.fimo.res$norm_score,first.fimo.res$pbm_id,max))
pbm[which(is.na(pbm))] <- 0
pbm <- matrix(data = pbm, nrow = length(unique(first.fimo.res$sequence_name)), ncol = length(unique(first.fimo.res$family)), byrow = TRUE,
              dimnames = list(sort(unique(first.fimo.res$sequence_name)),sort(unique(first.fimo.res$family))))
alltfs <- alltfs[which(alltfs$fam_id %in% colnames(pbm)),]
currlen <- length(unique(alltfs$tt_fam))
famkey <- famkey[which(famkey$fam_id %in% colnames(pbm)),]
dir.create('./modeling/', showWarnings = FALSE)
saveRDS(first.fimo.res,'./meme_stuff/fimo_res.rds')
saveRDS(alltfs, './modeling/tfkey.rds')
saveRDS(famkey, './modeling/famkey.rds')
saveRDS(pbm, "./modeling/pbm.rds")