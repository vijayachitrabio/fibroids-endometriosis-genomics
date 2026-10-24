suppressPackageStartupMessages({library(susieR);library(data.table)})
base <- "/Users/vijayachitramodhukur/Desktop/ECLAI/uterine_fibroids/outputs_v2/formal_ld_audit/susie_finemap"
out  <- "/Users/vijayachitramodhukur/Desktop/ECLAI/uterine_fibroids/Funciton_int_ge/revision_reruns/susie"
pips <- fread(file.path(base,"SuSiE_Shared_Loci_Variant_PIPs.csv")); summ <- fread(file.path(base,"SuSiE_Shared_Loci_Credible_Set_Summary.csv"))
cat("susieR", as.character(packageVersion("susieR")), "\n")
res <- list(); mem <- list(); diag <- list()
for (lc in unique(pips$locus)) for (tr in unique(pips[locus==lc]$trait)) {
  d <- pips[locus==lc & trait==tr]
  pre <- file.path(base,"ld_reference_subsets",paste0(lc,"_",tr))
  snps <- fread(paste0(pre,".snps.txt"),header=FALSE)$V1
  bim <- fread(paste0(pre,".ref.bim"),header=FALSE)
  R <- as.matrix(fread(paste0(pre,".ld.ld"),header=FALSE))
  if(nrow(R)!=nrow(bim)) {cat("dim mismatch",lc,tr,nrow(R),nrow(bim),"\n"); next}
  # align summary stats to the LD matrix order by SNP id; flip Z if alleles are swapped
  ix <- match(bim$V2, d$SNP); keep <- !is.na(ix)
  dd <- d[ix[keep]]; Rk <- R[keep,keep]; b <- bim[keep]
  flip <- ifelse(dd$EA==b$V6,1,ifelse(dd$EA==b$V5,-1,NA)); ok <- !is.na(flip)
  z <- dd$Z[ok]*flip[ok]; Rk <- Rk[ok,ok]; dd <- dd[ok]
  n <- summ[locus==lc & trait==tr]$n_eff[1]
  cat(lc,tr,"variants",length(z),"allele-mismatch dropped",sum(!ok),"n_eff",round(n),"\n")
  # concordance diagnostic: leave-one-out expected z (kriging) vs observed; lambda estimate
  lam <- tryCatch(estimate_s_rss(z,Rk,n=n),error=function(e) NA)
  kr <- tryCatch(kriging_rss(z,Rk,n=n),error=function(e) NULL)
  nout <- if(!is.null(kr)) sum(kr$conditional_dist$logLR>2 & abs(kr$conditional_dist$z_std_diff)>3,na.rm=TRUE) else NA
  diag[[length(diag)+1]] <- data.table(locus=lc,trait=tr,n_variants=length(z),alleles_dropped=sum(!ok),n_eff=n,estimated_s=lam,
      n_suspicious_variants_kriging=nout)
  if(!is.null(kr)) fwrite(data.table(SNP=dd$SNP,kr$conditional_dist)[abs(z_std_diff)>3 | logLR>2], file.path(out,paste0("kriging_flagged_",lc,"_",tr,".csv")))
  for (L in c(10,5,3,1)) {
    fit <- tryCatch(susie_rss(z,Rk,n=n,L=L,max_iter=500,coverage=0.95,min_abs_corr=0.5),error=function(e){cat("fail",lc,tr,L,conditionMessage(e),"\n");NULL})
    if(is.null(fit)) next
    cs <- fit$sets$cs
    res[[length(res)+1]] <- data.table(locus=lc,trait=tr,L=L,n_credible_sets=length(cs),converged=fit$converged,
        cs_sizes=paste(lengths(cs),collapse=","),lead_variants=if(length(cs)) paste(vapply(cs,function(i) dd$SNP[i[which.max(fit$pip[i])]],""),collapse=",") else "")
    if (L==10 && length(cs)) for (k in seq_along(cs)) mem[[length(mem)+1]] <- data.table(locus=lc,trait=tr,credible_set=names(cs)[k],
        SNP=dd$SNP[cs[[k]]],BP=dd$BP[cs[[k]]],EA=dd$EA[cs[[k]]],OA=dd$OA[cs[[k]]],marginal_P=dd$P[cs[[k]]],Z=z[cs[[k]]],PIP=fit$pip[cs[[k]]],
        cs_size=length(cs[[k]]),purity_min_abs_corr=fit$sets$purity$min.abs.corr[k])
  }
}
fwrite(rbindlist(res),file.path(out,"SuSiE_L_sensitivity_rerun.csv")); fwrite(rbindlist(mem),file.path(out,"SuSiE_credible_set_membership_L10.csv")); fwrite(rbindlist(diag),file.path(out,"SuSiE_LD_concordance_diagnostics.csv"))
cat("done\n")
