require(rstan)
require(R.utils)
rstan_options(auto_write = TRUE)

source("utils.R")

STANGLM=stan_model(file="betabinomial_glm.stan", save_dso=F, auto_write=F)

STANGLMNULL=stan_model(file="betabinomial_glm_null.stan", save_dso=F, auto_write=F)

# ys: numerator counts
# ns: denominator counts
# x: covariate vector
# Prior on concentration parameter is Gamma(concShape,concRate)
betaBinomialGLM_uni=function(ys,ns,x,concShape=1.0001,concRate=1e-4,...) {

  # try to get sensible initialization
  rat=ys/ns # ratios
  # moment estimator of concentration parameter
  conc=pmin( 1/var(rat, na.rm=T), 1000 ) 
  # thresholded moment estimator of the mean
  m=pmin( pmax( mean(rat, na.rm=T), 1/1000 ), 1-1/1000)
  init=list(conc=conc, mu=log(m/(1.0-m)))
    
  dat=list(N=length(ys),ys=ys,ns=ns,x=x,concShape=concShape,concRate=concRate)
  
  # fit the null model
  optimizing(STANGLMNULL, data=dat, init=init, algorithm="BFGS", hessian=F, as_vector=F, iter=1000, ...) -> stanresNull
  
  # initialize the alternative model using the fitted null model
  init=list(conc=stanresNull$par$conc, mu=stanresNull$par$mu, beta=0)
  optimizing(STANGLM, data=dat, init=init, algorithm="BFGS", hessian=T, iter=1000, ...) -> stanresFull
  
  #variance=tryCatch( {chol2inv(chol(-stanresFull$hessian))}, finally = { solve(-stanresFull$hessian) } )
  # get the Hessian at the mode
  variance=robustSolve(-stanresFull$hessian)
  dimnames(variance)=dimnames(stanresFull$hessian)
  # Laplace's estimate of the SE for beta
  betase=sqrt(variance["beta","beta"])
  # Wald z-score
  zscore=stanresFull$par["beta"]/betase
  # log likelihood ratio
  loglr=stanresFull$value - stanresNull$value
  
  list(beta=c(co=stanresFull$par["beta"],se=betase,p=2.0*pnorm(-abs(zscore))), loglr=loglr, lrtp=pchisq( 2.0*loglr, lower.tail = F , df=1 ), concUnderFull=stanresFull$par["conc"], concUnderNull=stanresNull$par$conc )
    
}

# numers: clusters x samples, numerator counts
# denom: cluster x samples, denominator counts
# x: covariate vector (assumed two class here), length samples
# minReads: minimum number of reads required in denominator
# minSampsPerGroup: smallest allowable number of samples in a group
# 
# Return values:
# coefs: # clusters x 3 matrix. Column1: effect size. Column2: standard error. Column3: p-value. 
# overdispersionUnderNull: estimated concentration parameters under null models
# overdispersionUnderFull: estimated concentration parameters under alternative models
# status: for each cluster, either one of success, not enough valid samples, or error thrown
runAllClusters=function(numers,denom,x,minReads=20,minSampsPerGroup=8,...){

  #require(doMC)
  #registerDoMC(if (parallel::detectCores()==4) 7 else parallel::detectCores())
  cat("Run all clusters\n")

  zz <- file( "/dev/null", open = "wt")
  sink(zz)
  sink(zz, type = "message")

  results=foreach(i=1:nrow(numers)) %dopar% {
    ys=numers[i,]
    ns=denom[i,]
    usable=ns>=minReads
    ta=table(x[usable])
    touse=ns>0
    ns=ns[touse]
    xh=x[touse]
    ys=ys[touse]
    if (length(ta)==2 & min(ta) >= minSampsPerGroup){
      if (sum(ys>0)>=5 & sum( (ys/ns)<1 )>=5 ) { # check not all 0
        tryCatch({
            res <- evalWithTimeout( { betaBinomialGLM_uni(ys,ns,xh,...) }, timeout=10, onTimeout="silent" ) 
            if (is.null(res)) "timeout" else res
        }, error=function(g) paste("error thrown:",g) 
                 )
      } else "almost all ys/ns is 0 or 1"
    } else "not enough valid samples"
  }
  sink(type="message")
  sink()

  cat("Done running clusters... process results\n")
  
  list( coefs=t(sapply(results,function(g) if (is.character(g)) rep(NA,3) else g$beta)), 
       overdispersionUnderNull=sapply(results,function(g) if (is.character(g))  NA else g$concUnderNull), 
       overdispersionUnderFull=sapply(results,function(g) if (is.character(g))  NA else g$concUnderFull),
       lrtp=sapply(results,function(g) if (is.character(g))  NA else g$lrtp[1]),
       status=sapply(results,function(g) if (is.character(g)) g else "success"))
}


testbbglm=function() {
  nsamp=100
  ns=rpois(nsamp,lambda=10)
  ys=rbinom(nsamp,ns,logistic(.3+.3*rnorm(nsamp)))
  betaBinomialGLM(ys,ns,runif(nsamp))
}

                                        # 
#optimizing(sm, data=list(N=N,ys=ys,ns=ns), init=init, algorithm="BFGS")
# 
# fit <- stan(model_code = schools_code, data = schools_dat, chains = 0)
# 
# n <- get_num_upars(fit)
# x <- runif(n, min = -2, max = 2)
# constrained <- constrain_pars(fit, x)
# is.list(constrained) # TRUE
# 
# my_pars <- relist(x^2, skeleton = constrained)
# unconstrained <- unconstrain_pars(fit, my_pars)
# log_prob(fit, unconstrained)
# grad_log_prob(fit, unconstrained)
