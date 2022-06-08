library(GenomicRanges)
library(rtracklayer)
add.info <- read.delim('../../mm39/info.txt', stringsAsFactors = FALSE)
mouse.gr.ucsc <- import("../../mm39/mm39.ncbiRefSeq.gtf")
for(n in 1:nrow(add.info)){
  end(mouse.gr.ucsc)[which(seqnames(mouse.gr.ucsc) == add.info[n,'Chrom'] & end(mouse.gr.ucsc) > add.info[n,'Pos'])] <- end(mouse.gr.ucsc)[which(seqnames(mouse.gr.ucsc) == add.info[n,'Chrom'] & end(mouse.gr.ucsc) > add.info[n,'Pos'])] + add.info[n,'Length']
  start(mouse.gr.ucsc)[which(seqnames(mouse.gr.ucsc) == add.info[n,'Chrom'] & start(mouse.gr.ucsc) > add.info[n,'Pos'])] <- start(mouse.gr.ucsc)[which(seqnames(mouse.gr.ucsc) == add.info[n,'Chrom'] & start(mouse.gr.ucsc) > add.info[n,'Pos'])] + add.info[n,'Length']
}
saveRDS(mouse.gr.ucsc,'mm39_gtf.rds')
export(mouse.gr.ucsc,'../../mm39/mm39_reformed.gtf', index = TRUE)
