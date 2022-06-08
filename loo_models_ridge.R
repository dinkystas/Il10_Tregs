library(ridge)
library(ggplot2)
library(ggrepel)
library(patchwork)
pbm <- readRDS('./modeling/pbm.rds')
atac.results <- readRDS('atac_results.rds')
numtfs <- ncol(pbm)
atac.results <- atac.results[rownames(pbm),]
svn.model <- readRDS('./modeling/svn_ridge_model.rds')
svp.model <- readRDS('./modeling/svp_ridge_model.rds')
svn.pvals <- readRDS('./modeling/svn_ridge_pvals.rds')
svp.pvals <- readRDS('./modeling/svp_ridge_pvals.rds')
TF.coeffs <- data.frame(TFs = names(coef(svn.model)[2:(1+numtfs)]), stringsAsFactors = FALSE)
TF.coeffs$svn_coefs <- coef(svn.model)[2:(1+numtfs)]
TF.coeffs$svp_coefs <- coef(svp.model)[2:(1+numtfs)]
TF.coeffs$svn_pvals <- svn.pvals$pval[,svn.pvals$chosen.nPCs]
TF.coeffs$svp_pvals <- svp.pvals$pval[,svp.pvals$chosen.nPCs]
sig.res <- subset(TF.coeffs, svn_pvals < 0.05 | svp_pvals < 0.05)
rownames(sig.res) <- sig.res$TFs
svp.models <- lapply(1:nrow(sig.res),function(i){
  currtf <- which(colnames(pbm) == sig.res$TFs[i])
  svpdata.mod <- data.frame(vals = atac.results[rownames(pbm),'lg2FC_svp'], pbm[,-currtf])
  model.mod <- linearRidge(vals ~ ., data = svpdata.mod)
})
svn.models <- lapply(1:nrow(sig.res),function(i){
  currtf <- which(colnames(pbm) == sig.res$TFs[i])
  svndata.mod <- data.frame(vals = atac.results[rownames(pbm),'lg2FC_svn'], pbm[,-currtf])
  model.mod <- linearRidge(vals ~ ., data = svndata.mod)
})
names(svp.models) <- sig.res$TFs
names(svn.models) <- sig.res$TFs
saveRDS(svp.models,'./modeling/svp_ridge_test_models.rds')
saveRDS(svn.models,'./modeling/svn_ridge_test_models.rds')