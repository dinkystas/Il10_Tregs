first.fimo.res <- readRDS('./meme_stuff/fimo_res.rds')
atac.results <- readRDS('atac_results.rds')
famkey <- readRDS('modeling/famkey.rds')
makebed <- function(family, filename){
  temp1 <- as.character(family)
  temp.motifs <- data.frame(motif = 1:length(which(first.fimo.res$family == temp1)),
                            peak = first.fimo.res$sequence_name[which(first.fimo.res$family == temp1)],
                            mstart = (first.fimo.res$start[which(first.fimo.res$family == temp1)] - 1),
                            mend = (first.fimo.res$stop[which(first.fimo.res$family == temp1)] - 1), stringsAsFactors = FALSE)
  temp.motifs <- merge(temp.motifs, atac.results[,c(1,2,3,5)], by = 'peak', sort = FALSE)
  temp.motifs$end <- temp.motifs$start + temp.motifs$mend
  temp.motifs$start <- temp.motifs$start + temp.motifs$mstart
  temp.motifs <- temp.motifs[,c(5,6,7,2)]
  temp.motifs$chr[which(temp.motifs$chr == 'chrX')] <- 'chr20'
  temp.motifs$chr[which(temp.motifs$chr == 'chrY')] <- 'chr21'
  temp.motifs$chr[which(temp.motifs$chr == 'chrM')] <- 'chr22'
  temp.motifs$chr <- unlist(lapply(1:nrow(temp.motifs),function(i){as.numeric(substring(temp.motifs$chr[i],4,5))}))
  temp.motifs <- temp.motifs[order(temp.motifs$chr, temp.motifs$start),]
  temp.motifs$chr <- paste0('chr',temp.motifs$chr)
  temp.motifs$chr[which(temp.motifs$chr == 'chr20')] <- 'chrX'
  temp.motifs$chr[which(temp.motifs$chr == 'chr21')] <- 'chrY'
  temp.motifs$chr[which(temp.motifs$chr == 'chr22')] <- 'chrM'
  rownames(temp.motifs) <- 1:nrow(temp.motifs)
  temp.motifs$motif <- paste(temp1,1:nrow(temp.motifs), sep = "_")
  write.table(temp.motifs, paste0('./Gen_out/Peak/', as.character(filename),'_motifs.bed'), sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
}


g1 <- unique(rna.results$Gene[which(rna.results$padj_svn < 0.05 & rna.results$RNA.cluster == 'Cluster_4')])
p1 <- names(which(pbm.tf[,'family_32']))
g2 <- unique(atac.results$Gene[which(atac.results$RNA.cluster == 'Cluster_4' & atac.results$peak %in% p1 & atac.results$lg2FC_svn > 0 & atac.results$padj_svn < 0.05)])

g3 <- unique(rna.results$Gene[which(rna.results$padj_svp < 0.05 & rna.results$RNA.cluster == 'Cluster_1')])
p2 <- names(which(pbm.tf[,'family_36']))
g4 <- unique(atac.results$Gene[which(atac.results$RNA.cluster == 'Cluster_1' & atac.results$peak %in% p2 & atac.results$lg2FC_svp < 0 & atac.results$padj_svp < 0.05)])
