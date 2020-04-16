#!/usr/bin/env Rscript

suppressMessages(library(optparse))
suppressMessages(library(DirichletMultinomial))
suppressMessages(library(TailRank))
suppressMessages(library(doMC))

#
# Load command line arguments and read in data
#

arguments <- parse_args(OptionParser(usage = "%prog [options] counts_file", description="LeafCutter outlier splicing command line tool. Required inputs:\n <counts_file>: Intron usage counts file. Must be .txt or .txt.gz, output from clustering pipeline.\n",option_list=list(
  make_option(c("-o","--output_prefix"), default = "leafcutter_outlier", help="The prefix for the output files,<prefix>_pVals.txt (containing p-value for each intron), <prefix>_clusterPvals.txt (containing p-value for each cluster) and <prefix>_effSize.txt (containing the effect sizes for each intron)  [default %default]"),
  make_option(c("-s","--max_cluster_size"), default=50, help="Don't test clusters with more introns than this [default %default]"), 
  make_option(c("-c","--min_coverage"), default=20, help="Require min_samples_per_group samples in each group to have at least this many reads [default %default]"), 
  make_option(c("-t","--timeout"), default=30, help="Maximum time (in seconds) allowed for a single optimization run [default %default]"),
  make_option(c("-p","--num_threads"), default=1, help="Number of threads to use [default %default]"))),
  positional_arguments = 1)

opt=arguments$opt
counts_file <- arguments$args[1]

cat("Loading counts from",counts_file,"\n")
if (!file.exists(counts_file)) stop("File ",counts_file," does not exist")
df <- read.table(counts_file, header=T, check.names = F)
introns <- rownames(df)
counts <- as.matrix(df)

if (opt$num_threads >1){
 registerDoMC(opt$num_threads)
}

cat("Settings:\n")
print(opt)


#
# Define all the functions for processing:
#

combineClusterPvals <- function(pVals){
    pVals <- as.numeric(pVals)
    return(pbeta(min(pVals),1,length(pVals),lower.tail = TRUE))
}

computeRightPvalBB <- function(count,totalCount,alpha,beta){
  #compute P(X>=count) = 1-P(X<=[count-1])
  if(count>0){
    val <- 1-pbb(count-1,totalCount,alpha,beta)
  }else{
    val <- 1
  }
  
  #fix numerical precision
  val <- max(.Machine$double.eps,min(val,1-.Machine$double.eps))
  return(val)
}


computeLeftPvalBB <- function(count,totalCount,alpha,beta){
  #compute P(X<=count) 
  val <- pbb(count,totalCount,alpha,beta)
  
  #fix numerical precision
  val <- max(.Machine$double.eps,min(val,1-.Machine$double.eps))
  return(val)
}

compute_effectSize_Pvalue <- function(count,totalCount,alphaEst,index,sumAlphaEst,clustMean){
  pRight <- computeRightPvalBB(count,totalCount,alphaEst[index],sumAlphaEst-alphaEst[index])
  pLeft <- computeLeftPvalBB(count,totalCount,alphaEst[index],sumAlphaEst-alphaEst[index])
  pVal <- min(1,2*min(pRight,pLeft))
  
  testVal <- (count/totalCount)
  effectSize <- testVal - clustMean[index]
  
  return(c(pVal=pVal,effectSize=effectSize,testVal=testVal))
}


runOutlierAnalysisCluster <- function(dataMat,alphaEst){
  #find cluster name
  cluID <- strsplit(colnames(dataMat)[1],":")[[1]][4]

  # initialize space
  effectSizesMat <- dataMat
  pValsMat <- dataMat
  testMat  <- dataMat
  
  # compute cluster summaries
  sumAlphaEst <- sum(alphaEst)
  clustMean <- alphaEst/sumAlphaEst
  
  # fill in matrices
  for(sample in seq_along(rownames(dataMat))){
    totalCount <- sum(dataMat[sample,])
    for(junction in seq_along(colnames(dataMat))){
      val <- compute_effectSize_Pvalue(dataMat[sample,junction],totalCount,alphaEst,junction,sumAlphaEst,clustMean)
      effectSizesMat[sample,junction] <- val["effectSize"]
      pValsMat[sample,junction] <- val["pVal"]
      testMat[sample,junction] <- val["testVal"]
    }
  }

  # compute combined pValue for the cluster
  combinedClustPval <- as.data.frame(apply(pValsMat,1,combineClusterPvals))
  colnames(combinedClustPval) <- cluID
  return(list(effectSizesMat=effectSizesMat,pValsMat=pValsMat,combinedClustPval=combinedClustPval,
              clustMean=clustMean,testMat=testMat))
}


#' Perform clusterwise splicing analysis.
#'
#' @param dataMat An [samples] x [introns] matrix of counts. The rownames must be of the form chr:start:end:cluid. If the counts file comes from the leafcutter clustering code this should be the case already.
#' @importFrom DirichletMultinomial dmn
#' @export
processCluster <- function(dataMat){
  #do laplace smoothing
  dataMat <- dataMat+1
  
  # model cluster
  dm_fit <- DirichletMultinomial::dmn(dataMat,1)
  
  # get results
  res <- runOutlierAnalysisCluster(dataMat,dm_fit@fit$Estimate)
}

writeResToFile <- function(results,output_prefix,verboseOutput=F){

  # print out first set with headers
  write.table( t(results[[1]]$pValsMat), 
               paste0(output_prefix,"_pVals.txt"), quote=F, col.names = T, row.names = T, sep="\t")
  write.table( t(results[[1]]$effectSizesMat), 
               paste0(output_prefix,"_effSize.txt"), quote=F, col.names = T, row.names = T, sep="\t")
  write.table( t(results[[1]]$combinedClustPval), 
               paste0(output_prefix,"_clusterPvals.txt"), quote=F, col.names = T, row.names = T, sep="\t")
  if(verboseOutput){
    write.table( results[[1]]$clustMean, 
                 paste0(output_prefix,"_controlMean.txt"), quote=F, col.names = F, row.names = T, sep="\t")
    write.table( t(results[[1]]$testMat), 
                 paste0(output_prefix,"_test.txt"), quote=F, col.names = T, row.names = T, sep="\t")
  }
  
  # append rest without headers
  for (index in setdiff(seq_along(results),1)){
    write.table( t(results[[index]]$pValsMat), 
                 paste0(output_prefix,"_pVals.txt"), quote=F, col.names = F, row.names = T, sep="\t", append = T)
    write.table( t(results[[index]]$effectSizesMat), 
                 paste0(output_prefix,"_effSize.txt"), quote=F, col.names = F, row.names = T, sep="\t", append = T)
    write.table( t(results[[index]]$combinedClustPval), 
                 paste0(output_prefix,"_clusterPvals.txt"), quote=F, col.names = F, row.names = T, sep="\t", append = T)
    if(verboseOutput){
      write.table( results[[index]]$clustMean, 
                   paste0(output_prefix,"_controlMean.txt"), quote=F, col.names = F, row.names = T, sep="\t", append = T)
      write.table( t(results[[index]]$testMat), 
                   paste0(output_prefix,"_test.txt"), quote=F, col.names = F, row.names = T, sep="\t", append = T)
    }
  }
}


#' Perform outlier splicing analysis.
#'
#' Parallelization across tested clusters is achieved using foreach/doMC, so the number of threads that will be used is determined by the cores argument passed to registerDoMC.
#'
#' @param counts An [introns] x [samples] matrix of counts. The rownames must be of the form chr:start:end:cluid. If the counts file comes from the leafcutter clustering code this should be the case already.
#' @param output_prefix The prefix for the output files.
#' @param max_cluster_size Don't test clusters with more introns than this
#' @param min_coverage Require min_samples_per_group samples in each group to have at least this many reads
#' @param timeout Maximum time (in seconds) allowed for a single optimization run
#' @return A per cluster list of results. Clusters that were not tested will be represented by a string saying why.
#' @import foreach
#' @importFrom R.utils withTimeout
#' @export
processAllClusters <- function(counts,output_prefix="leafcutter_outlier",max_cluster_size=10,min_coverage=20,timeout=10){
                              
  #Analyze cluster names
  splitIntrons <- do.call(rbind,strsplit(rownames(counts),":"))
  cluster_ids <- paste(splitIntrons[,1],splitIntrons[,4],sep = ":")
  cluster_sizes <- as.data.frame(table(cluster_ids))
  clu_names <- as.character(cluster_sizes$cluster_ids)
  cluster_sizes <- cluster_sizes$Freq
  names(cluster_sizes) <- clu_names
  
  results <- foreach (cluster_name=clu_names, .errorhandling = "remove") %dopar% {
    # filter out clusters
    if (cluster_sizes[cluster_name] > max_cluster_size)
      stop("Too many introns in cluster")
    junctions_in_cluster=which(cluster_ids==cluster_name)
    if (length(junctions_in_cluster) <= 1)
      stop("<=1 junction in cluster")
    cluster_counts <- t(counts[ junctions_in_cluster, ])
    sample_totals <- rowSums(cluster_counts)
    if (sum(sample_totals>=min_coverage)<=1)
      stop("<=1 sample with coverage>min_coverage")
    
    # process clusters that past filtering
    res <- R.utils::withTimeout( {
      processCluster(cluster_counts)
    }, timeout=timeout, onTimeout="error" )
    res$cluster_name <- cluster_name
    
    res
  } # end foreach
  
  writeResToFile(results,output_prefix)
  
  results
}



# Start Processing
cat("Running outlier splicing analysis...\n")

results <- processAllClusters(counts,
                              output_prefix=opt$output_prefix,
                              max_cluster_size=opt$max_cluster_size,
                              min_coverage=opt$min_coverage,
                              timeout=opt$timeout)

cat("Finished outlier splicing analysis.\n")



