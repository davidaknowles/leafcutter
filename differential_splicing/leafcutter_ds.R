#!/usr/bin/env Rscript
library(ptparse)
library(leafcutter)

arguments <- parse_args(OptionParser(usage = "%prog [options] counts_file groups_file", description="LeafCutter differential splicing command line tool. Required inputs:\n <counts_file>: Intron usage counts file. Must be .txt or .txt.gz, output from clustering pipeline.\n <groups_file>: Two column file: 1. sample names (must match column names in counts_file), 2. groups (currently only two groups, i.e. pairwise, supported. Some samples in counts_file can be missing from this file, in which case they will not be included in the analysis.",option_list=list(
  make_option(c("-o","--output_prefix"), default = "leafcutter_ds", help="The prefix for the two output files, <prefix>_cluster_significance.txt (containing test status, log likelihood ratio, degree of freedom, and p-value for each cluster) and <prefix>_effect_sizes.txt (containing the effect sizes for each intron)  [default %default]"),
  # make_option(c("-e","--confounders_file"), help="An optional file containing technical confounders (or principle components), where rows are the covariates and columns are the samples. The columns must match the columns in counts_file."),
  make_option(c("-s","--max_cluster_size"), default=Inf, help="Don't test clusters with more introns than this [default %default]"), 
  make_option(c("-i","--min_samples_per_intron"), default=5, help="Ignore introns used (i.e. at least one supporting read) in fewer than n samples [default %default]") , 
  make_option(c("-g","--min_samples_per_group"), default=3, help="Require this many samples in each group to have at least min_coverage reads [default %default]"), 
  make_option(c("-c","--min_coverage"), default=20, help="Require min_samples_per_group samples in each group to have at least this many reads [default %default]"), 
  make_option(c("-t","--timeout"), default=30, help="Maximum time (in seconds) allowed for a single optimization run [default %default]"),
  make_option(c("-p","--num_threads"), default=1, help="Number of threads to use [default %default]"))), 
  positional_arguments = 2)

opt=arguments$opt
counts_file=arguments$args[1]
groups_file=arguments$args[2]

cat("Loading counts from",counts_file,"\n")
if (!file.exists(counts_file)) stop("File ",counts_file," does not exist")
counts=read.table(counts_file, header=T)

cat("Loading metadata from",groups_file,"\n")
if (!file.exists(groups_file)) stop("File ",groups_file," does not exist")
meta=read.table(groups_file, header=F, stringsAsFactors = F)
colnames(meta)=c("sample","group")

counts=counts[,meta$sample]

meta$group=as.factor(meta$group)
group_names=levels(meta$group)

stopifnot(length(group_names)==2)

cat("Encoding as",group_names[1],"=0,",group_names[2],"=1\n")
numeric_x=as.numeric(meta$group)-1

require(doMC)
registerDoMC(opt$num_threads)

cat("Settings:\n")
print(opt)

cat("Running differential splicing analysis...\n")
results <- differential_splicing(counts, numeric_x, max_cluster_size=opt$max_cluster_size, min_samples_per_intron=opt$min_samples_per_intron, min_samples_per_group=opt$min_samples_per_group, min_coverage=opt$min_coverage, timeout=opt$timeout) 

cat("Saving results...\n")

write.table( cluster_results_table(results), paste0(opt$output_prefix,"_cluster_significance.txt"), quote=F, sep="\t", row.names = F)
write.table( leaf_cutter_effect_sizes(results), paste0(opt$output_prefix,"_effect_sizes.txt"), quote=F, col.names = F, sep="\t")

cat("All done, exiting\n")
