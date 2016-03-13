#' Perform splicing QTL analysis
#'
#' Parallelization across tested clusters is achieved using foreach/doMC, so the number of threads that will be used is determined by the cores argument passed to registerDoMC.
#'
#' @param counts An [introns] x [samples] matrix of counts. The rownames must be of the form chr:start:end:cluid. If the counts file comes from the leafcutter clustering code this should be the case already.
#' @param geno A [SNPs] x [samples] numeric matrix of the genotypes, typically encoded as 0,1,2, although in principle scaling shouldn't matter.
#' @param geno_meta SNP metadata, as a data.frame. Rows correspond to SNPs, must have a CHROM (with values e.g. chr15) and POS (position) column. 
#' @param pcs An optional [confounders] x [samples] matrix of technical confounders/PCs.
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
splicing_qtl=function(counts,geno,geno_meta,library_sizes,pcs=matrix(0,ncol(counts),0),snps_within=1e4,min_samples_per_intron=5,min_coverage=20,min_samples_per_group=8,debug=F,...) {
  
  introns=leafcutter:::get_intron_meta(rownames(counts))
  
  cluster_ids=paste(introns$chr,introns$clu,sep = ":")
  clusters_to_test=unique(cluster_ids)
  
  if (!debug)
     sink(file="/dev/null")
  res=foreach (clu=clusters_to_test, .errorhandling = if (debug) "stop" else "pass") %do% {
    
    print(clu)
    
    cluster_counts=t(counts[ cluster_ids==clu, ])
    sample_counts=rowSums(cluster_counts)
    samples_to_use=sample_counts>0
    
    if (sum(samples_to_use)<=1 | sum(sample_counts>=min_coverage)<=min_samples_per_group ) return("no samples_to_use")
    
    cluster_introns=introns[ cluster_ids %in% clu, ]
    m=mean(cluster_introns$middle)
    cis_snps = which( (abs( geno_meta$POS - m ) < snps_within) & (geno_meta$CHROM==cluster_introns$chr[1]) )
    
    introns_to_use=colSums(cluster_counts[samples_to_use,]>0)>=min_samples_per_intron
    cluster_counts=cluster_counts[,introns_to_use]

      sample_counts=sample_counts[samples_to_use]
    cluster_counts=cluster_counts[samples_to_use,]
    pcs_here=pcs[samples_to_use,,drop=F]

    cached_fit_null=NULL
    
    clures=foreach (cis_snp = cis_snps, .errorhandling = if (debug) "stop" else "pass") %dopar% {
      
      xh=as.numeric(geno[cis_snp,])
      
      if (length(unique(xh)) <= 1) return("Only one genotype")
      
      ta=table(xh[sample_counts>=min_coverage])

      if ( sum(ta >= min_samples_per_group) <= 1)
        return("not enough valid samples")

      if ( sum(introns_to_use)<2 )
        return("almost all ys/sample_counts is 0 or 1")
      
      linear_qtl(xh,cluster_counts,library_sizes)
    }
    
    setNames(clures,as.character(cis_snps))
  }
  if (!debug)
    sink()
  
  names(res)=clusters_to_test
  res 
}

linear_qtl=function(xh, counts, library_sizes) {
  # counts: columns = introns
  
  #counts=counts[ rowMeans(counts)>1, ]
  #counts=counts[ , colMeans(counts)>1 ]
  rat=sweep(as.matrix(counts),1,library_sizes,"/")
  
  x=data.frame(sample=rownames(counts), x=xh, stringsAsFactors = F)
  
  m=melt(rat)
  colnames(m)=c("sample", "intron", "rat")
  
  ij=inner_join(m, x, by="sample")
  
  #fit_null=MASS::glm.nb(count ~ sample + intron, data=ij)
  #fit_full= MASS::glm.nb(count ~ sample + intron + intron:x , data=ij)
  #anova(fit_null,fit_full)
  
  ij$log_count=log(ij$rat + .5)
  fit_null=lm(log_count ~ sample + intron, data=ij)
  fit_full=lm(log_count ~ sample + intron + intron:x , data=ij)
  anova(fit_null,fit_full)
}
