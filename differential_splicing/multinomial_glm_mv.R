
DIRICHLET_MULTINOMIAL_GLM=stan_model("dirichlet_multinomial_glm_mv_kM1.stan")

dirichlet_multinomial_glm <- function(x,y,concShape=1.0001,concRate=1e-4) {
  dat=list(N=nrow(x), K=ncol(y), P=ncol(x), y=y, x=x, concShape=concShape,concRate=concRate)
  stopifnot(nrow(x)==nrow(y))
  # fit model (optimization, don't judge)
  o=optimizing(DIRICHLET_MULTINOMIAL_GLM, data=dat, as_vector=F)

  list(value=o$value, conc=o$par$conc, beta=t(sweep( o$par$beta_raw - 1/K , 1, o$par$beta_scale, "*")))
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
