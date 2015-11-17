require(rstan)
# rstan_options(auto_write = TRUE)

source("utils.R")

STANGLM_MV=stan_model(file="betabinomial_glm_mv.stan", save_dso=F, auto_write=F)

STANGLM_FIX_CONC=stan_model(file="betabinomial_glm_fix_conc.stan", save_dso=F, auto_write=F)

# ys: numerator counts
# ns: denominator counts
# xFull: matrix of covariates for full model. First column must be ones. 
# xNull: matrix of covariates for null model. First column must be ones. 
# Prior on concentration parameter is Gamma(concShape,concRate)
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
  
  # TODO would like to use initialization but get error when length(beta)=1
  # Fit null model
  stanresNull <- optimizing(STANGLM_MV, data=dat, init=init, algorithm="BFGS", hessian=T, as_vector=F)
  
  # Initialize alternative model using null model
  betaInit=numeric(ncol(xFull))
  betaInit[1:ncol(xNull)]=stanresNull$par$beta
  initFull=list(conc=stanresNull$par$conc, beta=betaInit)
  
  # Fit alternative model
  datFull=dat
  datFull$x=xFull
  datFull$P=ncol(xFull)
  stanresFull <- optimizing(STANGLM_MV, data=datFull, init=initFull, algorithm="BFGS", hessian=T, as_vector=F)

  # Refit null model using concentration parameter from full fit
  # TODO: do we need this? 
  dat$conc=stanresFull$par$conc
  optimizing(STANGLM_FIX_CONC, data=dat, init=list(beta=stanresNull$par$beta), algorithm="BFGS" ) -> stanresNullFixConc
  
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

testbbglm=function() {
  nsamp=100
  ns=rpois(nsamp,lambda=10)
  ys=rbinom(nsamp,ns,logistic(.3+.3*rnorm(nsamp)))
  xNull=cbind(numeric(nsamp)+1)
  xFull=cbind(xNull,runif(nsamp))
  betaBinomialGLM(ys,ns,xFull,xNull)
}
