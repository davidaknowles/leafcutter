
#' Benjamini-Hochberg FDR correction
#' @param p P-values
#' @return Adjusted p-values (q-values)
#' @export
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

#' @import graphics
pqplot=function(p, ...) { pqplotHelper(p); abline(0,1) }

#' Plot multiple qq plot for p-values
#'
#' @param pvalues A named list of numeric vectors of p-values.
#' @return a ggplot
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

#' Make a data.frame of meta data about the introns
#' @param introns Names of the introns
#' @return Data.frame with chr, start, end, cluster id and "middle" 
#' @export
get_intron_meta=function(introns){
  intron_meta=do.call(rbind,strsplit(introns,":"))
  colnames(intron_meta)=c("chr","start","end","clu")
  intron_meta=as.data.frame(intron_meta,stringsAsFactors = F)
  intron_meta$start=as.numeric(intron_meta$start)
  intron_meta$end=as.numeric(intron_meta$end)
  intron_meta$middle=.5*(intron_meta$start+intron_meta$end) # TODO: can remove this now? 
  intron_meta
}

mahalanobis_outlier=function(x) {
    ei=eigen(cov(x))
    if (any(ei$values<=0.0)) return(numeric(nrow(x))+1)
    prec=ei$vectors %*% diag(1/ei$values) %*% t(ei$vectors)
    mah_dist=mahalanobis( x, colMeans(x), prec, inverted=T )
    pchisq(mah_dist, df=ncol(x), lower.tail=F) * nrow(x)
}

sanitize_simplex=function(x, eps=1e-6) {
  x[x<eps]=eps
  x[x>(1.0-eps)]=(1.0-eps)
  x/sum(x)
}


#' Work out which gene each cluster belongs to. Note the chromosome names used in the two inputs must match. 
#' @param intron_meta Data frame describing the introns, usually from get_intron_meta
#' @param exons_table Table of exons, see e.g. /data/gencode19_exons.txt.gz
#' @return Data.frame with cluster ids and genes separated by commas
#' @import dplyr
#' @export
map_clusters_to_genes=function(intron_meta, exons_table) {
    gene_df=foreach (chr=sort(unique(intron_meta$chr)), .combine=rbind) %dopar% {
    
        intron_chr=intron_meta[ intron_meta$chr==chr, ]
        exons_chr=exons_table[exons_table$chr==chr, ]

        exons_chr$temp=exons_chr$start
        intron_chr$temp=intron_chr$end
        three_prime_matches=inner_join( intron_chr, exons_chr, by="temp")

        exons_chr$temp=exons_chr$end
        intron_chr$temp=intron_chr$start
        five_prime_matches=inner_join( intron_chr, exons_chr, by="temp")

        all_matches=rbind(three_prime_matches, five_prime_matches)[ , c("clu", "gene_name")]

        all_matches=all_matches[!duplicated(all_matches),]
        
        if (nrow(all_matches)==0) return(NULL)
        all_matches$clu=paste(chr,all_matches$clu,sep=':')
        all_matches
    }

    clu_df=gene_df %>% group_by(clu) %>% summarize(genes=paste(gene_name, collapse = ","))
    
    class(clu_df)="data.frame"

    clu_df
}

#' Adding "chr" if it's not present
#' @param chrs Chromosome or cluster names
#' @return Data.frame with cluster ids and genes separated by commas
#' @export
add_chr=function(chrs)
    if (!grepl("chr",chrs[1])) paste0("chr",chrs) else chrs