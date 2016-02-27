library(leafcutter)

test_dirichlet_multinomial_glm_mc <- function() {

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
  
  an=dirichlet_multinomial_anova_mc(x,x[,1:2],y,robust=F)

  an=dirichlet_multinomial_anova_mc(x,x[,1:2],y,robust=T)
  
  # check result looks reasonable
  #plot( as.numeric(beta), as.numeric(dm$beta) )
  #abline(0,1)
}

test_dirichlet_multinomial_glm_mc()

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
  #plot( as.numeric(beta), as.numeric(dm$beta) )
  #abline(0,1)
}


testbbglm=function() {
  nsamp=100
  ns=rpois(nsamp,lambda=10)
  ys=rbinom(nsamp,ns,logistic(.3+.3*rnorm(nsamp)))
  xNull=cbind(numeric(nsamp)+1)
  xFull=cbind(xNull,runif(nsamp))
  betaBinomialGLM(ys,ns,xFull,xNull)
}

compare_bb_and_dm=function() {
  set.seed(1)
  nsamp=100
  ns=rpois(nsamp,lambda=10)
  ys=rbinom(nsamp,ns,logistic(.3+.3*rnorm(nsamp)))
  y=cbind(ys,ns-ys)
  xNull=cbind(numeric(nsamp)+1)
  xFull=cbind(xNull,runif(nsamp))
  
  dm_fit=dirichlet_multinomial_anova(xFull,xNull,y)
  bb_fit=betaBinomialGLM(ys,ns,xFull,xNull)
  
  stopifnot( abs( (dm_fit$fit_full$par$conc - bb_fit$concUnderFull)/dm_fit$fit_full$par$conc ) < 1e-2 )
  stopifnot( abs( (dm_fit$fit_null$par$conc - bb_fit$concUnderNull)/dm_fit$fit_null$par$conc ) < 1e-2 )
  stopifnot( abs( dm_fit$loglr - bb_fit$loglr[1] ) < 1e-3 )
}

