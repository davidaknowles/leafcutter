require(doMC)
require(R.utils)
registerDoMC(7)
source("~/Dropbox/splicing/leafcutter/differential_splicing/multinomial_glm_multi_conc.R",echo=T)

cluster_results_table=function(results) {
  rows=foreach(res=results) %do% 
  { if ( !is.list(res) ) NULL else 
    c(loglr=res$loglr, df=res$df, p=res$lrtp) }
  names(rows)=names(results)
  as.data.frame(do.call(rbind, rows))
}

leafcutter_status=function(results) {
  unlist( foreach(res=results) %do% { if ( !is.list(res) ) res else "Successfully tested" } )
}

beta_real=function(r) 
  sweep(r$beta_raw - 1.0/ncol(r$beta_raw), 1, r$beta_scale, "*") 

leaf_cutter_effect_sizes=function(results) {
  unlist( foreach(res=results) %do% {
    if (is.list(res)) beta_real( res$fit_full$par )[2,] else NULL
  } )
}

differential_splicing=function(counts, x, max_cluster_size=10, min_samples_per_intron=5, min_samples_per_group=4, min_coverage=20, timeout=10) {
  
  intron_names=rownames(counts)
  introns=as.data.frame(do.call(rbind,strsplit(intron_names,":")))
  colnames(introns)=c("chr","start","end","clu")
  cluster_ids=paste(introns$chr,introns$clu,sep = ":")
  
  cluster_sizes=as.data.frame(table(cluster_ids))
  clu_names=as.character(cluster_sizes$cluster_ids)
  cluster_sizes=cluster_sizes$Freq
  names(cluster_sizes)=clu_names

  zz <- file( "/dev/null", open = "wt")
  sink(zz)
  sink(zz, type = "message")
  
  results=foreach (cluster_name=clu_names) %dopar% {
    
  #results=lapply(clu_names, function(cluster_name) {
    if (cluster_sizes[cluster_name] > max_cluster_size)
      return("Too many introns in cluster")
    cluster_counts=t(counts[ cluster_ids==cluster_name, ])
    sample_totals=rowSums(cluster_counts)
    samples_to_use=sample_totals>0 
    if (sum(samples_to_use)<=1 ) 
      return("<=1 sample with coverage>0")
    sample_totals=sample_totals[samples_to_use]
    if (sum(sample_totals>=min_coverage)<=1) 
      return("<=1 sample with coverage>min_coverage")
    x_subset=x[samples_to_use]
    cluster_counts=cluster_counts[samples_to_use,]
    introns_to_use=colSums(cluster_counts>0)>=min_samples_per_intron # only look at introns used by at least 5 samples
    if ( sum(introns_to_use)<2 ) 
      return("<2 introns used in >=min_samples_per_intron samples")
    cluster_counts=cluster_counts[,introns_to_use]
    ta=table(x_subset[sample_totals>=min_coverage])
    if (sum(ta >= min_samples_per_group)<2) # at least two groups have this (TODO: continuous x)
      return("Not enough valid samples") 
    xFull=cbind(1,x_subset)
    xNull=xFull[,1,drop=F]
    #tryCatch({
      res <- evalWithTimeout( { 
        dirichlet_multinomial_anova_mc(xFull,xNull,cluster_counts)
      }, timeout=timeout, onTimeout="silent" ) 
      if (is.null(res)) "timeout" else res
    #}, error=function(g) g)
  }
  
  sink(type="message")
  sink()
  
  statuses=leafcutter_status(results)
  
  cat("Differential splicing summary:\n")
  print(as.data.frame(table(statuses)))
  
  names(results)=clu_names
  results
}



