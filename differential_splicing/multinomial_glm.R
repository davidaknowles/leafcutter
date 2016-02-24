require(rstan)

DIRICHLET_MULTINOMIAL_GLM=stan_model("dirichlet_multinomial_glm_mv_kM1.stan")

dirichlet_multinomial_glm <- function(x,y,concShape=1.0001,concRate=1e-4) {
  dat=list(N=nrow(x), K=ncol(y), P=ncol(x), y=y, x=x, concShape=concShape,concRate=concRate)
  stopifnot(nrow(x)==nrow(y))

  o=optimizing(DIRICHLET_MULTINOMIAL_GLM, data=dat, as_vector=F)
  list(value=o$value, conc=o$par$conc, beta=t(sweep( o$par$beta_raw - 1/dat$K , 1, o$par$beta_scale, "*")))
}

dirichlet_multinomial_anova <- function(xFull,xNull,y,concShape=1.0001,concRate=1e-4) {
  K=ncol(y)
  dat=list(N=nrow(xNull), K=K, P=ncol(xNull), y=y, x=xNull, concShape=concShape,concRate=concRate)
  # fit model (optimization, don't judge)
  fit_null=optimizing(DIRICHLET_MULTINOMIAL_GLM, data=dat, as_vector=F)
  
  dat_full=list(N=nrow(xFull), K=K, P=ncol(xFull), y=y, x=xFull, concShape=concShape,concRate=concRate)
  # beta_raw is PxK
  init=list(beta_raw=matrix(1/K,ncol(xFull),K), beta_scale=rep(1,ncol(xFull)), conc=fit_null$par$conc)
  init$beta_raw[1:ncol(xNull),]=fit_null$par$beta_raw
  init$beta_scale[1:ncol(xNull)]=fit_null$par$beta_scale
  # fit model (optimization, don't judge)
  fit_full=optimizing(DIRICHLET_MULTINOMIAL_GLM, data=dat_full, init=init, as_vector=F)
  
  loglr=fit_full$value-fit_null$value
  df=( ncol(xFull)-ncol(xNull) )*(K-1)
  list( loglr=loglr, df=df, lrtp=pchisq( 2.0*loglr, lower.tail = F , df=df ), fit_null=fit_null, fit_full=fit_full)
}

test_dirichlet_multinomial_glm <- function() {

  N=100 # number of samples
  depth=50 # number in each sample
  K=4 # number of classes
  P=3 # number of covariates

  softmax=function(g) exp(g)/sum(exp(g))

  beta=matrix(rnorm(K*P),K,P) # true coefficient matrix
  x=lapply(1:N, function(g) runif(P)) # covariates
  # observed counts
  y=lapply(x, function(g) rmultinom(1, depth, softmax(beta %*% g))) 
  y=t(do.call(cbind,y)) 
  x=do.call(rbind,x)
  dm=dirichlet_multinomial_glm(x,y)
  
  # check result looks reasonable
  plot( as.numeric(beta), as.numeric(dm$beta) )
  abline(0,1)
}
