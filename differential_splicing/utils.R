require(ggplot2)

bh=function(p) { q=p; q[!is.na(p)]=p.adjust(pmin(p[!is.na(p)],1),method="fdr"); q }
bfc=function(p) min(p,na.rm=T)*sum(!is.na(p))

# logistic (sigmoid) function
logistic=function(g) 1/(1+exp(-g))

# get inverse of PSD matrix using eigendecomposition
# make robust by thresholding all eigenvalues to be > eigenThreshold
robustSolve=function(si, eigenThreshold=0.01){
  svdSigma=eigen(si)
  svdSigma$values=pmax(svdSigma$values,eigenThreshold)
  svdSigma$vectors %*% diag(1/svdSigma$values) %*% t(svdSigma$vectors)
}

fisherCombined=function(p) pchisq( -2*sum(log(p),na.rm=T), df=2*sum(!is.na(p)), lower.tail = F)

pqplotHelper=function(p, ...) qqplot(-log10(runif(length(p))), -log10(p), pch=16, ...)

pqplot=function(p, ...) { pqplotHelper(p); abline(0,1) }

multiqq=function(bfcp) {
  punif=-log10(runif(max(sapply(bfcp,length))))
  df=do.call(rbind, foreach( i=seq_len(length(bfcp))) %do%  {
    df=as.data.frame( qqplot(punif[1:length(bfcp[[i]])], -log10(bfcp[[i]]), plot.it=F) )
    df$npcs=names(bfcp)[i]
    df
  } ) 
  df$npcs=factor(df$npcs, names(bfcp))
  ggplot(df, aes(x,y,col=npcs)) + geom_point() + geom_abline(intercept=0,slope=1) + theme_bw(base_size=18) + xlab("Expected -log10(p)") + ylab("Observed -log10(p)") 
}