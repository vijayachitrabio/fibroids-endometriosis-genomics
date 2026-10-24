suppressPackageStartupMessages({library(data.table);library(LAVA)})
a <- commandArgs(TRUE); info<-a[1]; ov<-a[2]; tag<-a[3]; i0<-as.integer(a[4]); i1<-as.integer(a[5])
ref <- "/Users/vijayachitramodhukur/Desktop/ECLAI/uterine_fibroids/mr_analysis_2026_04_16/reference_data/g1000_eur"
t0<-Sys.time(); input <- process.input(input.info.file=info, sample.overlap.file=if(ov=="none") NULL else ov, ref.prefix=ref, input.dir=getwd())
cat("input s:",as.numeric(Sys.time()-t0,units="secs"),"\n"); loci <- read.loci("lava_loci_2495_blocks.tsv"); loci <- loci[i0:min(i1,nrow(loci)),]
out <- list(); qc <- list()
for(i in seq_len(nrow(loci))){
  loc <- loci[i,]
  locus <- tryCatch(process.locus(loc,input,drop.failed=FALSE,max.block.size=3000),error=function(e) NULL)
  if(is.null(locus)){qc[[i]]<-data.table(LOC=loc$LOC,status="process_locus_null");next}
  res <- tryCatch(run.univ.bivar(locus,univ.thresh=1,p.values=TRUE,CIs=TRUE),error=function(e) NULL)
  if(is.null(res)||is.null(res$biv)||nrow(res$biv)==0){
    u <- if(!is.null(res$univ)) as.data.table(res$univ) else data.table()
    qc[[i]]<-data.table(LOC=loc$LOC,status="no_bivariate",n_snps=locus$n.snps,K=locus$K,
      univ_p=if(nrow(u)) paste(signif(as.numeric(u$p),3),collapse=";") else NA_character_,h2=if(nrow(u)) paste(signif(as.numeric(u$h2.obs),3),collapse=";") else NA_character_); next}
  u <- as.data.table(res$univ)
  b <- as.data.table(res$biv)
  out[[i]] <- cbind(data.table(LOC=loc$LOC,CHR=loc$CHR,START=loc$START,STOP=loc$STOP,n_snps=locus$n.snps,K=locus$K,
     univ_p_endo=u$p[1],univ_p_uf=u$p[2],h2_endo=u$h2.obs[1],h2_uf=u$h2.obs[2]),b)
  qc[[i]]<-data.table(LOC=loc$LOC,status="ok",n_snps=locus$n.snps,K=locus$K)
}
fwrite(rbindlist(out,fill=TRUE),paste0("res_",tag,".csv")); fwrite(rbindlist(qc,fill=TRUE),paste0("qc_",tag,".csv"))
cat("done",tag,"\n")
