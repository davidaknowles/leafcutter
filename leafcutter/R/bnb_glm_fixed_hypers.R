#model_to_use=stan_model("~/Dropbox/splicing/leafcutter/leafcutter/exec/bnb_glm_reparam.stan")

#' Beta negative binomial GLM likelihood ratio test for a single cluster
#'
#' @param x [samples] vector
#' @param y [samples] x [introns] matrix of intron usage counts
#' @param concShape Gamma shape parameter for concentration parameter
#' @param concShape Gamma rate parameter for concentration parameter
#' @param fit_null Optionally the fitted null model (used in \code{\link{splicing_qtl}} to save repeatedly fitting the null for each cis-SNP)
#' @importFrom rstan optimizing
#' @importFrom reshape2 melt
#' @export
bnb_glm_fixed_hypers <- function(x, cluster_counts, hypers, fit_null=NULL, ...) {
  
  model_to_use=leafcutter:::stanmodels$bnb_glm_fixed_hypers
  
  melted_counts=reshape2::melt(data.frame(cluster_counts, sample=rownames(cluster_counts), x=x, check.names = F), id.vars=c("sample","x"), value.name = "count")
  colnames(melted_counts)[3]="intron"
  
  rownames(hypers)=colnames(cluster_counts)
  hypers_n=hypers[as.character(melted_counts$intron),]
  
  xNull=model.matrix( ~ sample + intron - 1, data=melted_counts)
  #det( t(xNull) %*% xNull )
  
  dat_null=list(N=nrow(xNull), P=ncol(xNull), y=melted_counts$count, x=xNull, r=hypers_n$r, b=hypers_n$b )
  #dat_null=list(N=nrow(xNull), P=ncol(xNull), y=melted_counts$count, x=xNull, sqrtrbShape=1.01,sqrtrbRate=0.001, sqrtrOverBShape=3, sqrtrOverBRate=3 )
  if (is.null(fit_null)) fit_null=rstan::optimizing(model_to_use, data=dat_null, as_vector=F, ...)
  
  intron_mm=model.matrix( ~ intron , data=melted_counts )
  intron_mm=intron_mm[,2:ncol(intron_mm),drop=F] # remove intercept
  
  x_new=sweep(intron_mm,1,x,"*")
  
  xFull=cbind( xNull, x_new )
  #det( t(xFull) %*% xFull )
  
  dat_full=dat_null
  dat_full$P=ncol(xFull)
  dat_full$x=xFull
  
  df=dat_full$P - dat_null$P
  
  init=fit_null$par
  init$beta=c(init$beta, numeric(df))

  fit_full=rstan::optimizing(model_to_use, data=dat_full, init=init, as_vector=F, ...)
  
  loglr=fit_full$value-fit_null$value
  
  pchisq( 2.0*loglr, lower.tail = F , df=df )
  
  refit_null_flag=F
  
  lrtp=pchisq( 2.0*loglr, lower.tail = F , df=df )
  if (lrtp < .001) {
    init=fit_full$par
    init$beta=init$beta[seq_len(dat_null$P)]

    refit_null=rstan::optimizing(model_to_use, data=dat_null, init=init, as_vector=F, ...)
    if (refit_null$value > fit_null$value) {
      refit_null_flag=T
      fit_null=refit_null
      loglr=fit_full$value-fit_null$value
    }
  }
  
  list( loglr=loglr, df=df, lrtp=pchisq( 2.0*loglr, lower.tail = F , df=df ), fit_null=fit_null, fit_full=fit_full, refit_null_flag=refit_null_flag)
}
