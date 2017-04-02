
#'Dirichlet multinomial GLM (single overdispersion parameter)
#'
#'We recommend using \code{\link{dirichlet_multinomial_anova_mc}} instead. 
#' 
#' @param x [samples] x [covariates] matrix
#' @param y [samples] x [introns] matrix of intron usage counts
#' @param concShape Gamma shape parameter for concentration parameter
#' @param concRate Gamma rate parameter for concentration parameter
#' @importFrom rstan optimizing
#' @export
dirichlet_multinomial_glm <- function(x,y,concShape=1.0001,concRate=1e-4) {
  dat=list(N=nrow(x), K=ncol(y), P=ncol(x), y=y, x=x, concShape=concShape,concRate=concRate)
  stopifnot(nrow(x)==nrow(y))

  o=rstan::optimizing(stanmodels$dm_glm, data=dat, as_vector=F)
  list(value=o$value, conc=o$par$conc, beta=t(sweep( o$par$beta_raw - 1/dat$K , 1, o$par$beta_scale, "*")))
}

#' Dirichlet multinomial GLM likelihood ratio test for a single cluster (single overdispersion parameter)
#'
#'We recommend using \code{\link{dirichlet_multinomial_glm_anova_mc}} instead. 
#'
#' @param xFull [samples] x [covariates] matrix for the alternative model
#' @param xNull [samples] x [covariates] matrix for the null model
#' @param y [samples] x [introns] matrix of intron usage counts
#' @param concShape Gamma shape parameter for concentration parameter
#' @param concRate Gamma rate parameter for concentration parameter
#' @importFrom rstan optimizing
#' @export
dirichlet_multinomial_anova <- function(xFull,xNull,y,concShape=1.0001,concRate=1e-4) {
  K=ncol(y)
  dat=list(N=nrow(xNull), K=K, P=ncol(xNull), y=y, x=xNull, concShape=concShape,concRate=concRate)

  fit_null=rstan::optimizing(stanmodels$dm_glm, data=dat, as_vector=F)
  
  dat_full=list(N=nrow(xFull), K=K, P=ncol(xFull), y=y, x=xFull, concShape=concShape,concRate=concRate)
  # beta_raw is PxK
  init=list(beta_raw=matrix(1/K,ncol(xFull),K), beta_scale=rep(1,ncol(xFull)), conc=fit_null$par$conc)
  init$beta_raw[1:ncol(xNull),]=fit_null$par$beta_raw
  init$beta_scale[1:ncol(xNull)]=fit_null$par$beta_scale

  fit_full=rstan::optimizing(stanmodels$dm_glm, data=dat_full, init=init, as_vector=F)
  
  loglr=fit_full$value-fit_null$value
  df=( ncol(xFull)-ncol(xNull) )*(K-1)
  list( loglr=loglr, df=df, lrtp=pchisq( 2.0*loglr, lower.tail = F , df=df ), fit_null=fit_null, fit_full=fit_full)
}

