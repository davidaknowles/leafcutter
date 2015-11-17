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
runAllClusters=function(numers,x,minReads=20,minSampsPerGroup=8,...){
  
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
        status=sapply(results,function(g) if (is.character(g)) g else "success"))
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
