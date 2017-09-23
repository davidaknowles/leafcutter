

#' Dirichlet multinomial GLM for confounder removal and PSI quantification
#'
#' @param x [samples] x [covariates] matrix of confounders
#' @param cluster_counts [samples] x [introns] matrix of intron usage counts
#' @param concShape Gamma shape parameter for concentration parameter
#' @param concRate Gamma rate parameter for concentration parameter
#' @param debug Whether to give verbose output from rstan.
#' @param init Can be one of {"smart", "random"}. smart uses an method of moments estimator to get a reasonable initialization. The seed for "random" can be set through the ... arguments passed to rstan::optimizing.
#' @param smart_init_regularizer Used to protect against colinear covariates. 
#' @param residual_sigma Std of regularization on residuals.
#' @param ... will be passed on the rstan::optimizing, so can be used for example to set the algorithm used (default is LBFGS) or the random seed if random initialization is requested. 
#' @importFrom rstan optimizing
#' @import foreach
#' @import dplyr
#' @import magrittr
#' @export
quantify_psi_one_cluster <- function(cluster_counts,x,protected,concShape=1.0001,concRate=1e-4, debug=F, init_strategy="smart", smart_init_regularizer=0.001, residual_sigma=10, ...) {
  K=ncol(cluster_counts)
  
  model_to_use=stanmodels$dm_glm_multi_conc
  
  dat_null=list(N=nrow(x), K=K, P=ncol(x), y=cluster_counts, x=x, concShape=concShape,concRate=concRate)
  
  if (init_strategy=="smart") {
    y_norm = sweep(cluster_counts+1,1,rowSums(cluster_counts+1),"/") %>% log()
    beta_mm = solve( t(x) %*% x + smart_init_regularizer * diag(ncol(x)), t(x) %*% y_norm )
    beta_norm = sweep(beta_mm, 1, rowMeans(beta_mm), "-")
    beta_scale = foreach(i=seq_len(nrow(beta_norm)), .combine = c) %do% {
      beta_row = beta_norm[i,]
      up=beta_row[which.max(abs(beta_row))] / (1-1/K)
      down=beta_row[which.max(-sign(up)*beta_row)] / (1/K)
      up-down
    }
    beta_raw=sweep(beta_norm,1,beta_scale+1e-20,"/") + 1/K
    beta_raw=sweep(beta_raw,1,rowSums(beta_raw),"/")
    init=list(beta_scale=array(beta_scale), beta_raw=beta_raw, conc=rep(10.0,K))
  }
  
  fit=optimizing(model_to_use, data=dat_null, init=init, as_vector=F, verbose=debug, ...)

  dimnames(fit$par$beta_raw)=list(colnames(x),colnames(cluster_counts))
  
  dat_null$beta=t(beta_real(fit$par))
  dat_null$conc=fit$par$conc
  dat_null$residual_sigma=residual_sigma
  
  fit_psi=optimizing(stanmodels$dm_glm_mc_psi, data=dat_null, init=0, as_vector=F, verbose=debug, ...)
  
  normalize=function(g) { g/sum(g) }
  softmax=function(g) normalize(exp(g))
  
  lo_residuals = fit_psi$par$residual + x[,protected,drop=F] %*% t(dat_null$beta[,protected,drop=F])
  
  lo_residuals %>%
    apply( 1,softmax) %>% # implicit t() here
    sweep(1,fit$par$conc,"*") %>% 
    apply( 2,normalize) %>% 
    set_rownames(colnames(cluster_counts)) %>%
    set_colnames(rownames(cluster_counts))
}

#' Confounder removal and PSI quantification
#'
#' @param counts An [introns] x [samples] matrix of counts. The rownames must be of the form chr:start:end:cluid. If the counts file comes from the leafcutter clustering code this should be the case already.
#' @param x A [samples] x [confounders] numeric matrix to be controlled for in the GLM. Factors should already have been converted to a 1-of-(K-1) encoding, e.g. using model.matrix (see scripts/leafcutter_ds.R for how to do this). 
#' @param protected Indices (boolean or integer) corresponding to which columns of x should be protected rather than regressed out. 
#' @param timeout Maximum time (in seconds) allowed for a single optimization run
#' @param debug If true writes more output
#' @param init One of 'smart' (default) or 'random'. If 'random' you can pass an additional arg "seed" for reproducibility. 
#' @return A per cluster list of results. Clusters that were not tested will be represented by a string saying why.
#' @import foreach
#' @importFrom R.utils evalWithTimeout
#' @export
quantify_psi=function(counts, x, protected, timeout=10, debug=F, init="smart", ...) {
  
  introns=get_intron_meta(rownames(counts))
  cluster_ids=paste(introns$chr,introns$clu,sep = ":")
  
  cluster_sizes=as.data.frame(table(cluster_ids))
  clu_names=as.character(cluster_sizes$cluster_ids)
  
  if (!debug) {
    zz <- file( "/dev/null", open = "wt")
    sink(zz)
    sink(zz, type = "message")
  }
  
  tryCatch( {
  psi_matrix=do.call(rbind, foreach (cluster_name=clu_names, .errorhandling = if (debug) "stop" else "pass") %dopar% {
    cluster_counts=t(counts[ cluster_ids==cluster_name, ])
    evalWithTimeout( { 
      quantify_psi_one_cluster(cluster_counts, x, protected, debug=debug, init=init,...)
    }, timeout=timeout, onTimeout=if (debug) "silent" else "warning" ) 
  } )
  }, error=function(g) {
    if (!debug) {
      sink(type="message")
      sink()
    }
    print(g)
    return(NULL)
  } )
  
  if (!debug) {
    sink(type="message")
    sink()
  }
  
  psi_matrix
}

