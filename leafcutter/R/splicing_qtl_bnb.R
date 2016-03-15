#' Perform splicing QTL analysis
#'
#' Parallelization across tested clusters is achieved using foreach/doMC, so the number of threads that will be used is determined by the cores argument passed to registerDoMC.
#'
#' @param counts An [introns] x [samples] matrix of counts. The rownames must be of the form chr:start:end:cluid. If the counts file comes from the leafcutter clustering code this should be the case already.
#' @param geno A [SNPs] x [samples] numeric matrix of the genotypes, typically encoded as 0,1,2, although in principle scaling shouldn't matter.
#' @param geno_meta SNP metadata, as a data.frame. Rows correspond to SNPs, must have a CHROM (with values e.g. chr15) and POS (position) column. 
#' @param snps_within Window from center of cluster in which to test SNPs. 
#' @param max_cluster_size Don't test clusters with more introns than this
#' @param min_samples_per_intron Ignore introns used (i.e. at least one supporting read) in fewer than n samples
#' @param min_samples_per_group Require this many samples in each group to have at least min_coverage reads
#' @param min_coverage Require min_samples_per_group samples in each group to have at least this many reads
#' @param timeout Maximum time (in seconds) allowed for a single optimization run
#' @param debug Turn on to see output from rstan.
#' @return A per cluster list of results. For each cluster this is a list over tested SNPs. SNPs that were not tested will be represented by a string saying why.
#' @import foreach
#' @importFrom R.utils evalWithTimeout
#' @export
splicing_qtl_bnb=function(counts,geno,geno_meta,snps_within=1e4,min_samples_per_intron=5,min_coverage=20,min_samples_per_group=8,timeout=10,debug=F,...) {
  
  introns=leafcutter:::get_intron_meta(rownames(counts))
  
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
    #m=mean(cluster_introns$middle)
    #cis_snps = which( (abs( geno_meta$POS - m ) < snps_within) & (geno_meta$CHROM==cluster_introns$chr[1]) )
    cis_snps = which(  ( (min(cluster_introns$start) - snps_within) <  geno_meta$POS ) & ( geno_meta$POS < (max(cluster_introns$end) + snps_within)) & (geno_meta$CHROM==cluster_introns$chr[1]) )
    
    introns_to_use=colSums(cluster_counts[samples_to_use,]>0)>=min_samples_per_intron
    cluster_counts=cluster_counts[,introns_to_use]

    if (sum(introns_to_use)<=1) return("<=1 usable introns")

      sample_counts=sample_counts[samples_to_use]
    cluster_counts=cluster_counts[samples_to_use,]
    #pcs_here=pcs[samples_to_use,,drop=F]

    cached_fit_null=NULL
    
    clures=foreach (cis_snp = cis_snps, .errorhandling = if (debug) "stop" else "pass") %do% {
      
      xh=as.numeric(geno[cis_snp,])
      
      if (length(unique(xh)) <= 1) return("Only one genotype")
      
      ta=table(xh[sample_counts>=min_coverage])

      if ( sum(ta >= min_samples_per_group) <= 1)
        return("not enough valid samples")

      if ( sum(introns_to_use)<2 )
        return("almost all ys/sample_counts is 0 or 1")

        if (debug & !is.null(cached_fit_null)) cat("Using cached null fit.\n")
        res <- R.utils::evalWithTimeout( { bnb_glm(xh,cluster_counts,fit_null=cached_fit_null,...) }, timeout=timeout, onTimeout="silent" )
        if (is.null(res)) "timeout" else {
          cached_fit_null=res$fit_null
          res
        }
    }
    
    names(clures)=as.character(cis_snps)
    clures
  }
  if (!debug)
    sink()
  
  names(res)=clusters_to_test
  res 
}
