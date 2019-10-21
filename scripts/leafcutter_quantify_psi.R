#!/usr/bin/env Rscript
library(optparse)

arguments <- parse_args(OptionParser(usage = "%prog [options] counts_file", description="LeafCutter PSI quantification command line tool. Required inputs:\n <counts_file>: Intron usage counts file. Must be .txt or .txt.gz, output from clustering pipeline.\n <confounders_file>: ",option_list=list(
  make_option(c("-c","--confounders_file"), default=NULL, help="One+K column file: 1. sample names (must match column names in counts_file), 2. Additional columns are used to specify confounders, e.g. batch/sex/age. Numeric columns will be treated as continuous, so use e.g. batch1, batch2, batch3 rather than 1, 2, 3 if you want a categorical variable."),
  make_option(c("-o","--output_file"), default = "leafcutter_psi.txt.gz", help="The output file, will be gzipped [default %default]"),
   make_option(c("-t","--timeout"), default=30, help="Maximum time (in seconds) allowed for a single optimization run [default %default]"),
  make_option(c("-p","--num_threads"), default=1, help="Number of threads to use [default %default]"),
  make_option(c("--init"), default="smart", help="One of 'smart' (default) or 'random'."), 
  make_option(c("--seed"), default=12345, help="Random seed if using random initialization."))),
  positional_arguments = 1)

library(leafcutter)

opt=arguments$opt
counts_file=arguments$args[1]

cat("Loading counts from",counts_file,"\n")
if (!file.exists(counts_file)) stop("File ",counts_file," does not exist")
counts=read.table(counts_file, header=T, check.names = F)

if (!is.null(opt$confounders_file)) {
  cat("Loading counfounders from",opt$confounders_file,"\n")
  if (!file.exists(opt$confounders_file)) stop("File ",opt$confounders_file," does not exist")
  meta=read.table(opt$confounders_file, header=F, stringsAsFactors = F)
  colnames(meta)[1]="sample"
  counts=counts[,meta$sample]
  
  confounders=meta[,2:ncol(meta),drop=F]
  # scale continuous confounders
  for (i in seq_len(ncol(confounders)))
    if (is.numeric(confounders[,i]))
      confounders[,i]=scale(confounders[,i])
  # convert factors to one-of-K encoding
  confounders=model.matrix( ~., data=confounders )
  confounders=confounders[,2:ncol(confounders),drop=F] # remove intercept
  
} else {
  confounders = matrix(0,nrow=ncol(counts),ncol=0)
}

require(doMC)
registerDoMC(opt$num_threads)

cat("Settings:\n")
print(opt)

cat("Running PSI quantification\n")
x=cbind(intercept=1,confounders)
psi_matrix <- quantify_psi(counts, x, protected=1, timeout=opt$timeout, init=opt$init, seed=opt$seed , debug = FALSE) 

stopifnot(!is.null(psi_matrix))

cat("Saving results...\n")
gz=gzfile(opt$output_file, "wb")
write.table(psi_matrix, gz, quote=F, sep="\t",col.names = T,row.names = T)
close(gz)

cat("All done, exiting\n")
