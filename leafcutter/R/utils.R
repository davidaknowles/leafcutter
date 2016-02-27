
bh=function(p) { q=p; q[!is.na(p)]=p.adjust(pmin(p[!is.na(p)],1),method="fdr"); q }

bfc=function(p) min(p,na.rm=T)*sum(!is.na(p))

#' logistic (sigmoid) function
#' @param g The log odds
#' @return 1/(1+e^-g)
#' @export
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

#' Plot multiple qq plot for p-values
#'
#' @param pvales A list of numeric vectors of p-values.
#' @return ggplot
#' @import ggplot2
#' @export
multiqq=function(pvalues) {
  punif=-log10(runif(max(sapply(pvalues,length))))
  df=do.call(rbind, foreach( i=seq_len(length(pvalues))) %do%  {
    df=as.data.frame( qqplot(punif[1:length(pvalues[[i]])], -log10(pvalues[[i]]), plot.it=F) )
    df$group=names(pvalues)[i]
    df
  } ) 
  df$group=factor(df$group, names(pvalues))
  ggplot(df, aes(x,y,col=group)) + geom_point() + geom_abline(intercept=0,slope=1) + theme_bw(base_size=18) + xlab("Expected -log10(p)") + ylab("Observed -log10(p)") 
}

get_intron_meta=function(introns){
  intron_meta=do.call(rbind,strsplit(introns,":"))
  colnames(intron_meta)=c("chr","start","end","clu")
  intron_meta=as.data.frame(intron_meta,stringsAsFactors = F)
  intron_meta$start=as.numeric(intron_meta$start)
  intron_meta$end=as.numeric(intron_meta$end)
  intron_meta$middle=.5*(intron_meta$start+intron_meta$end)
  intron_meta
}

mahalanobis_outlier=function(x) {
    ei=eigen(cov(x))
    if (any(ei$values<=0.0)) return(numeric(nrow(x))+1)
    prec=ei$vectors %*% diag(1/ei$values) %*% t(ei$vectors)
    mah_dist=mahalanobis( x, colMeans(x), prec, inverted=T )
    pchisq(mah_dist, df=ncol(x), lower.tail=F) * nrow(x)
}
