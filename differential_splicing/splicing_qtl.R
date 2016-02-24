source("dm_glm_multi_conc.R")

get_intron_meta=function(introns){
  intron_meta=do.call(rbind,strsplit(introns,":"))
  colnames(intron_meta)=c("chr","start","end","clu")
  intron_meta=as.data.frame(intron_meta,stringsAsFactors = F)
  intron_meta$start=as.numeric(intron_meta$start)
  intron_meta$end=as.numeric(intron_meta$end)
  intron_meta$middle=.5*(intron_meta$start+intron_meta$end)
  intron_meta
}

splicing_qtl=function(counts,geno,geno_meta,pcs=matrix(0,ncol(counts),0),permute=F,snps_within=1e4,min_samples_per_intron=5,min_coverage=20,min_samples_per_group=8,timeout=10,debug=F) {
  
  introns=get_intron_meta(rownames(counts))
  
  cluster_ids=paste(introns$chr,introns$clu,sep = ":")
  clusters_to_test=unique(cluster_ids)
  
  if (!debug)
     sink(file="/dev/null")
  res=foreach (clu=clusters_to_test, .errorhandling = if (debug) "stop" else "pass") %dopar% {
    
    print(clu)
    
    cluster_counts=t(counts[ cluster_ids==clu, ])
    sample_counts=rowSums(cluster_counts)
    samples_to_use=sample_counts>0
    
    if (sum(samples_to_use)<=1 | sum(sample_counts>=min_coverage)<=min_samples_per_group ) return("no samples_to_use")
    
    cluster_introns=introns[ cluster_ids %in% clu, ]
    m=mean(cluster_introns$middle)
    cis_snps = which( abs( geno_meta$POS - m ) < snps_within )
    
    sample_counts=sample_counts[samples_to_use]
    
    cluster_counts=cluster_counts[samples_to_use,]
    introns_to_use=colSums(cluster_counts>0)>=min_samples_per_intron
    cluster_counts=cluster_counts[,introns_to_use]
    
    pcs_here=pcs[samples_to_use,,drop=F]
    
    cached_fit_null=NULL
    
    clures=foreach (cis_snp = cis_snps, .errorhandling = if (debug) "stop" else "pass") %do% {
      
      xh=as.numeric(geno[cis_snp,samples_to_use])
      
      if (length(unique(xh)) <= 1) return("Only one genotype")
      
      if (permute) {
        possible_perms=foreach(temptemp=1:1000) %do% sample(xh)
        perm_correlations=foreach(pp=possible_perms, .combine=c) %do% cor(xh,pp)
        xh=possible_perms[[which.min(abs(perm_correlations))]]
      }
      
      ta=table(xh[sample_counts>=min_coverage])

      if ( sum(ta >= min_samples_per_group) <= 1)
        return("not enough valid samples")

      if ( sum(introns_to_use)<2 )
        return("almost all ys/sample_counts is 0 or 1")
      
      #tryCatch({
        xFull=cbind(1,pcs_here,xh)
        xNull=cbind(1,pcs_here)
        if (debug & !is.null(cached_fit_null)) cat("Using cached null fit.\n")
        res <- evalWithTimeout( { dirichlet_multinomial_anova_mc(xFull,xNull,cluster_counts,fit_null=cached_fit_null) }, timeout=timeout, onTimeout="silent" )
        if (is.null(res)) "timeout" else {
          cached_fit_null=res$fit_null
          res
        }
      #}, error=function(g) as.character(g) )
    }
    
    names(clures)=as.character(cis_snps)
    clures
  }
  if (!debug)
    sink()
  
  names(res)=clusters_to_test
  res 
}
