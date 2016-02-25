
#' Dirichlet multinomial GLM
#'
#' @param x [samples] x [covariates] matrix
#' @param y [samples] x [introns] matrix of intron usage counts
#' @param concShape Gamma shape parameter for concentration parameter
#' @param concShape Gamma rate parameter for concentration parameter
#' @importFrom rstan optimizing
#' @export
dirichlet_multinomial_glm_mc <- function(x,y,concShape=1.0001,concRate=1e-4) {
  dat=list(N=nrow(x), K=ncol(y), P=ncol(x), y=y, x=x, concShape=concShape,concRate=concRate)
  stopifnot(nrow(x)==nrow(y))

  o=rstan::optimizing(stanmodels$dm_glm_multi_conc, data=dat, as_vector=F)
  list(value=o$value, conc=o$par$conc, beta=t(sweep( o$par$beta_raw - 1/dat$K , 1, o$par$beta_scale, "*")))
}

#' Dirichlet multinomial GLM likelihood ratio test for a single cluster
#'
#' @param xFull [samples] x [covariates] matrix for the alternative model
#' @param xNull [samples] x [covariates] matrix for the null model
#' @param y [samples] x [introns] matrix of intron usage counts
#' @param concShape Gamma shape parameter for concentration parameter
#' @param concShape Gamma rate parameter for concentration parameter
#' @param fit_null Optionally the fitted null model (used in \code{\link{splicing_qtl}} to save repeatedly fitting the null for each cis-SNP)
#' @importFrom rstan optimizing
#' @export
dirichlet_multinomial_anova_mc <- function(xFull,xNull,y,concShape=1.0001,concRate=1e-4, fit_null=NULL) {
  K=ncol(y)
  
  dat_null=list(N=nrow(xNull), K=K, P=ncol(xNull), y=y, x=xNull, concShape=concShape,concRate=concRate)
  
  # fit null model
  if (is.null(fit_null)) fit_null=rstan::optimizing(stanmodels$dm_glm_multi_conc, data=dat_null, as_vector=F)

  colnames(fit_null$par$beta_raw)=colnames(y)
  rownames(fit_null$par$beta_raw)=colnames(xNull)
  
  dat_full=list(N=nrow(xFull), K=K, P=ncol(xFull), y=y, x=xFull, concShape=concShape,concRate=concRate)
  # beta_raw is PxK
  init=list(beta_raw=matrix(1e-4,ncol(xFull),K), beta_scale=rep(1,ncol(xFull)), conc=fit_null$par$conc)

  # beta_raw must live _in_ the simplex
  beta_raw_sanitized=fit_null$par$beta_raw
  beta_raw_sanitized[beta_raw_sanitized<1e-6]=1e-6
  beta_raw_sanitized[beta_raw_sanitized>(1.0-1e-6)]=(1.0-1e-6)
  init$beta_raw[1:ncol(xNull),]=beta_raw_sanitized
  init$beta_raw=sweep(init$beta_raw, 1, rowSums(init$beta_raw), "/") 
  
  init$beta_scale[1:ncol(xNull)]=fit_null$par$beta_scale
  stopifnot(all(is.finite(unlist(init))))
  # fit fit model
  fit_full=rstan::optimizing(stanmodels$dm_glm_multi_conc, data=dat_full, init=init, as_vector=F)

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
    refit_null=rstan::optimizing(stanmodels$dm_glm_multi_conc, data=dat_null, init=init, as_vector=F)
    if (refit_null$value > fit_null$value) {
      refit_null_flag=T
      fit_null=refit_null
      loglr=fit_full$value-fit_null$value
    }
  }
  
  list( loglr=loglr, df=df, lrtp=pchisq( 2.0*loglr, lower.tail = F , df=df ), fit_null=fit_null, fit_full=fit_full, refit_null_flag=refit_null_flag)
}
