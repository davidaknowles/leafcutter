#' Beta binomial GLM likelihood ratio testing
#'
#' We recommend using \code{\link{dirichlet_multinomial_glm_anova}} instead, but this is here for comparison.
#' 
#' @param ys numerator counts
#' @param ns denominator counts
#' @param xFull matrix of covariates for full model. First column must be ones. 
#' @param xNull matrix of covariates for null model. First column must be ones. 
#' @param concShape Gamma shape parameter for concentration parameter
#' @param concRate Gamma rate parameter for concentration parameter
#' @param ... will be passed on the rstan::optimizing, so can be used for example to set the algorithm used (default is LBFGS).
#' @importFrom rstan optimizing
#' @export
#' @import stats
betaBinomialGLM=function(ys,ns,xFull,xNull,concShape=1.0001,concRate=1e-4,...) {
  stopifnot(all(xNull==xFull[,1:ncol(xNull)]))
  stopifnot(all(xNull[,1]==1))
  # try to get sensible initialization
  rat=ys/ns # ratios
  # moment estimator of concentration parameter
  conc=pmin( 1/var(rat, na.rm=T), 1000 ) 
  # thresholded moment estimator of the mean
  m=pmin( pmax( mean(rat, na.rm=T), 1/1000 ), 1-1/1000)
  
  betaInit=numeric(ncol(xNull))
  betaInit[1]=log(m/(1.0-m))
  init=list(conc=conc, beta=as.array(betaInit))
  #sampling(STANGLM, data=list(N=length(ys),P=ncol(xNull),ys=ys,ns=ns,x=xNull), chains=1) -> stanresNull
  
  dat=list(N=length(ys),P=ncol(xNull),ys=ys,ns=ns,x=xNull,concShape=concShape,concRate=concRate)
  
  # Fit null model
  stanresNull <- rstan::optimizing(stanmodels$bb_glm, data=dat, init=init, algorithm="BFGS", hessian=T, as_vector=F)
  
  # Initialize alternative model using null model
  betaInit=numeric(ncol(xFull))
  betaInit[1:ncol(xNull)]=stanresNull$par$beta
  initFull=list(conc=stanresNull$par$conc, beta=betaInit)
  
  # Fit alternative model
  datFull=dat
  datFull$x=xFull
  datFull$P=ncol(xFull)
  stanresFull <- rstan::optimizing(stanmodels$bb_glm, data=datFull, init=initFull, algorithm="BFGS", hessian=T, as_vector=F)

  # Refit null model using concentration parameter from full fit
  # TODO: do we need this? 
  dat$conc=stanresFull$par$conc
  rstan::optimizing(stanmodels$bb_glm_fix_conc, data=dat, init=list(beta=stanresNull$par$beta), algorithm="BFGS" ) -> stanresNullFixConc
  
  # Function to extract coefficients, standard errors, and Wald p-values
  getCoefs=function(res,x) {
      variance=robustSolve(-res$hessian)
      dimnames(variance)=dimnames(res$hessian)
      betase=sqrt(diag(variance))[paste0("beta.",1:ncol(x))]
      beta=res$par[paste0("beta[",1:ncol(x),"]")]
      zscore=res$par$beta / betase
      data.frame(co=res$par$beta,se=betase,p=2.0*pnorm(-abs(zscore)))
  }
  
  # Two approaches to calculating the log likelihood statistic
  # 1. Null+alternative model are allowed different concentration parameters
  # 2. Both use concentration parameter from alternative (full) model.
  loglr=c( stanresFull$value - stanresNull$value, 
      stanresFull$value - stanresNullFixConc$value )
  
  list(betaNull=getCoefs(stanresNull,xNull), betaFull=getCoefs(stanresFull,xFull), loglr=loglr, lrtp=pchisq( 2.0*loglr, lower.tail = F , df=ncol(xFull)-ncol(xNull) ), concUnderFull=stanresFull$par$conc, concUnderNull=stanresNull$par$conc )
}

