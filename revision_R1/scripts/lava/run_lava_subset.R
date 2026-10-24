# Pass 2: run LAVA on a given block list, writing each block as it finishes (checkpointed, resumable)
suppressPackageStartupMessages({library(data.table);library(LAVA)})
a <- commandArgs(TRUE); info<-a[1]; ov<-a[2]; tag<-a[3]; blocks<-a[4]
ref <- "/Users/vijayachitramodhukur/Desktop/ECLAI/uterine_fibroids/mr_analysis_2026_04_16/reference_data/g1000_eur"
input <- process.input(input.info.file=info, sample.overlap.file=if(ov=="none") NULL else ov, ref.prefix=ref, input.dir=getwd())
loci <- read.loci("lava_loci_2495_blocks.tsv"); loci <- loci[loci$LOC %in% readLines(blocks),]
rf <- paste0("res_",tag,".csv"); qf <- paste0("qc_",tag,".csv")
done <- if(file.exists(qf)) fread(qf)$LOC else character()
for(i in seq_len(nrow(loci))){
  loc <- loci[i,]; if(loc$LOC %in% done) next
  locus <- tryCatch(process.locus(loc,input,drop.failed=FALSE,max.block.size=3000),error=function(e) e)
  if(inherits(locus,"error")||is.null(locus)){fwrite(data.table(LOC=loc$LOC,status="process_locus_failed",msg=if(inherits(locus,"error")) conditionMessage(locus) else "NULL"),qf,append=file.exists(qf));next}
  res <- tryCatch(run.univ.bivar(locus,univ.thresh=1,p.values=TRUE,CIs=TRUE),error=function(e) e)
  if(inherits(res,"error")){fwrite(data.table(LOC=loc$LOC,status="univ_bivar_error",msg=conditionMessage(res)),qf,append=file.exists(qf));next}
  u <- as.data.table(res$univ)
  if(is.null(res$biv)||nrow(res$biv)==0){fwrite(data.table(LOC=loc$LOC,status="no_bivariate",msg=paste(signif(as.numeric(u$p),3),collapse=";")),qf,append=file.exists(qf));next}
  b <- as.data.table(res$biv)
  fwrite(cbind(data.table(LOC=loc$LOC,CHR=loc$CHR,START=loc$START,STOP=loc$STOP,n_snps=locus$n.snps,K=locus$K,
     univ_p_endo=u$p[1],univ_p_uf=u$p[2],h2_endo=u$h2.obs[1],h2_uf=u$h2.obs[2]),b),rf,append=file.exists(rf))
  fwrite(data.table(LOC=loc$LOC,status="ok",msg=""),qf,append=file.exists(qf))
}
cat("done",tag,"\n")
