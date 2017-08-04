#' Parse output of differential_splicing
#'
#' Parse output of \code{\link{differential_splicing}} and make a per cluster results table
#'
#' @param results From \code{\link{differential_splicing}}
#' @return Data.frame with columns status, log likelihood ratio, degrees of freedom, p-value
#' @export
cluster_results_table=function(results) {
  as.data.frame(cbind(cluster=names(results), foreach(res=results, .combine=rbind) %do% 
{ if ( is.character(res) | ("error" %in% class(res)) ) data.frame(status=as.character(res), loglr=NA, df=NA, p=NA) else 
  data.frame(status="Success", loglr=res$loglr, df=res$df, p=res$lrtp) } ) )
}

#' Convert K-1 representation of parameters to real
#'
#' @param r A parameter object from fitting dm_glm or dm_glm_multi_conc. Must have members beta_raw and beta_scale. 
#' @return Coefficient matrix.
#' @export
beta_real=function(r) 
  sweep(r$beta_raw - 1.0/ncol(r$beta_raw), 1, r$beta_scale, "*") 

#' Per intron effect sizes
#'
#' @param results From \code{\link{differential_splicing}}
#' @return data.frame with a row for every tested intron and columns: intron, log effect size, baseline proportions, proportions in the second condition, and resulting deltaPSI
#' @importFrom dplyr bind_rows
#' @export
leaf_cutter_effect_sizes=function(results) {
  softmax=function(g) exp(g)/sum(exp(g))
  all_introns=foreach(res=results, .combine = bind_rows) %do% {
    if ( is.character(res) | ("error" %in% class(res)) ) NULL else {
       beta=beta_real( res$fit_full$par )
       data.frame( intron=colnames(beta), logef=beta[2,], baseline=softmax(beta[1,]), perturbed=softmax(beta[1,]+beta[2,]), stringsAsFactors = F )
    }      
  } 
  all_introns$deltapsi=with(all_introns, perturbed-baseline)
  all_introns
}


#' Perform pairwise differential splicing analysis.
#'
#' Parallelization across tested clusters is achieved using foreach/doMC, so the number of threads that will be used is determined by the cores argument passed to registerDoMC.
#'
#' @param counts An [introns] x [samples] matrix of counts. The rownames must be of the form chr:start:end:cluid. If the counts file comes from the leafcutter clustering code this should be the case already.
#' @param x A [samples] numeric vector, should typically be 0s and 1s, although in principle scaling shouldn't matter.
#' @param confounders A [samples] x [confounders] numeric matrix to be controlled for in the GLM. Factors should already have been converted to a 1-of-(K-1) encoding, e.g. using model.matrix (see scripts/leafcutter_ds.R for how to do this). Can be NULL, implying no covariates are controlled for. 
#' @param max_cluster_size Don't test clusters with more introns than this
#' @param min_samples_per_intron Ignore introns used (i.e. at least one supporting read) in fewer than n samples
#' @param min_samples_per_group Require this many samples in each group to have at least min_coverage reads
#' @param min_coverage Require min_samples_per_group samples in each group to have at least this many reads
#' @param timeout Maximum time (in seconds) allowed for a single optimization run
#' @param robust Whether to use the robust model (explicitly models outliers). Generally not required/recommended for differential splicing.
#' @param debug If true writes more output
#' @return A per cluster list of results. Clusters that were not tested will be represented by a string saying why.
#' @import foreach
#' @importFrom R.utils evalWithTimeout
#' @export
differential_splicing=function(counts, x, confounders=NULL, max_cluster_size=10, min_samples_per_intron=5, min_samples_per_group=4, min_coverage=20, timeout=10, robust=F, debug=F) {
  
  stopifnot(ncol(counts)==length(x))
  
  introns=get_intron_meta(rownames(counts))
  cluster_ids=paste(introns$chr,introns$clu,sep = ":")
  
  cluster_sizes=as.data.frame(table(cluster_ids))
  clu_names=as.character(cluster_sizes$cluster_ids)
  cluster_sizes=cluster_sizes$Freq
  names(cluster_sizes)=clu_names
  
  if (!debug) {
    zz <- file( "/dev/null", open = "wt")
    sink(zz)
    sink(zz, type = "message")
  }
  
  results=foreach (cluster_name=clu_names, .errorhandling = "pass") %dopar% {
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
    if (!is.null(confounders)) {
        ch=confounders[samples_to_use,,drop=F]
        ch=ch[ , apply(ch,2,sd)>0.0, drop=F ]
        xFull=cbind(xFull,ch)
        xNull=cbind(xNull,ch)
    }
    res <- R.utils::evalWithTimeout( { 
      dirichlet_multinomial_anova_mc(xFull,xNull,cluster_counts,robust=robust, debug=debug)
    }, timeout=timeout, onTimeout="silent" ) 
    if (is.null(res)) "timeout" else res
  }
  
  if (!debug) {
    sink(type="message")
    sink()
  }
  
    names(results)=clu_names

    statuses=cluster_results_table(results)$status
  
  cat("Differential splicing summary:\n")
  print(as.data.frame(table(statuses)))
  
  results
}

