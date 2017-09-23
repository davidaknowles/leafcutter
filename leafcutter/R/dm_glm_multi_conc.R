

#' Dirichlet multinomial GLM likelihood ratio test for a single cluster
#'
#' @param xFull [samples] x [covariates] matrix for the alternative model
#' @param xNull [samples] x [covariates] matrix for the null model
#' @param y [samples] x [introns] matrix of intron usage counts
#' @param concShape Gamma shape parameter for concentration parameter
#' @param concRate Gamma rate parameter for concentration parameter
#' @param robust Whether to include an outlier model (i.e. use dm_glm_multi_conc_robust rather than dm_glm_multi_conc)
#' @param outlier_prior_a Only used for robust model. The outlier probability outlier_prob ~ Beta(outlier_prior_a,outlier_prior_b)
#' @param outlier_prior_b Only used for robust model. The outlier probability outlier_prob ~ Beta(outlier_prior_a,outlier_prior_b)
#' @param fit_null Optionally cache the fitted null model to save repeatedly fitting the null for each cis-SNP when sQTL mapping)
#' @param debug Whether to give verbose output from rstan.
#' @param init Can be one of {"smart", "random"}. smart uses an method of moments estimator to get a reasonable initialization. The seed for "random" can be set through the ... arguments passed to rstan::optimizing.
#' @param smart_init_regularizer Used to protect against colinear covariates. 
#' @param ... will be passed on the rstan::optimizing, so can be used for example to set the algorithm used (default is LBFGS) or the random seed if random initialization is requested. 
#' @importFrom rstan optimizing
#' @import foreach
#' @import dplyr
#' @export
dirichlet_multinomial_anova_mc <- function(xFull,xNull,y,concShape=1.0001,concRate=1e-4, robust=T, outlier_prior_a=1.01, outlier_prior_b=100, fit_null=NULL, debug=F, init="smart", smart_init_regularizer=0.001, ...) {
  K=ncol(y)
  
  model_to_use=if (robust) stanmodels$dm_glm_robust else stanmodels$dm_glm_multi_conc
  
  dat_null=list(N=nrow(xNull), K=K, P=ncol(xNull), y=y, x=xNull, concShape=concShape,concRate=concRate, outlier_prior_a=outlier_prior_a, outlier_prior_b=outlier_prior_b)
  
  if (init=="smart") {
    y_norm = sweep(y+1,1,rowSums(y+1),"/") %>% log()
    beta_mm = solve( t(xNull) %*% xNull + smart_init_regularizer * diag(ncol(xNull)), t(xNull) %*% y_norm )
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
  
  # fit null model
  if (is.null(fit_null)) 
    fit_null=rstan::optimizing(model_to_use, data=dat_null, init=init, as_vector=F, verbose=debug, ...)

  colnames(fit_null$par$beta_raw)=colnames(y)
  rownames(fit_null$par$beta_raw)=colnames(xNull)
  
  dat_full=dat_null
  dat_full$P=ncol(xFull)
  dat_full$x=xFull

  # beta_raw is PxK
  init=list(beta_raw=matrix(1e-4,ncol(xFull),K), beta_scale=rep(1,ncol(xFull)), conc=fit_null$par$conc, outlier_prob=fit_null$par$outlier_prob)

  # beta_raw must live _in_ the simplex
  beta_raw_sanitized=fit_null$par$beta_raw
  beta_raw_sanitized[beta_raw_sanitized<1e-6]=1e-6
  beta_raw_sanitized[beta_raw_sanitized>(1.0-1e-6)]=(1.0-1e-6)
  init$beta_raw[1:ncol(xNull),]=beta_raw_sanitized
  init$beta_raw=sweep(init$beta_raw, 1, rowSums(init$beta_raw), "/") 
  
  init$beta_scale[1:ncol(xNull)]=fit_null$par$beta_scale
  stopifnot(all(is.finite(unlist(init))))
  # fit fit model
  fit_full=rstan::optimizing(model_to_use, data=dat_full, init=init, as_vector=F, verbose=debug,...)

  colnames(fit_full$par$beta_raw)=colnames(y)
  rownames(fit_full$par$beta_raw)=colnames(xFull)
  
  loglr=fit_full$value-fit_null$value
  df=( ncol(xFull)-ncol(xNull) )*(K-1)
  
  refit_null_flag=F
  
  lrtp=pchisq( 2.0*loglr, lower.tail = F , df=df )
  if (lrtp < .001) {
    init=fit_full$par
    init$beta_raw=init$beta_raw[seq_len(dat_null$P),,drop=F]
    init$beta_raw[init$beta_raw<1e-6]=1e-6
    init$beta_raw[init$beta_raw>(1.0-1e-6)]=(1.0-1e-6)
    init$beta_raw=sweep(init$beta_raw, 1, rowSums(init$beta_raw), "/") 
    init$beta_scale=as.array(init$beta_scale[seq_len(dat_null$P)])
    refit_null=rstan::optimizing(model_to_use, data=dat_null, init=init, as_vector=F, verbose=debug,...)
    if (refit_null$value > fit_null$value) {
      refit_null_flag=T
      fit_null=refit_null
      loglr=fit_full$value-fit_null$value
    }
  }
  
  list( loglr=loglr, df=df, lrtp=pchisq( 2.0*loglr, lower.tail = F , df=df ), fit_null=fit_null, fit_full=fit_full, refit_null_flag=refit_null_flag)
}
