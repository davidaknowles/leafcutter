#!/usr/bin/env Rscript
library(optparse)
library(leafcutter)

arguments <- parse_args(OptionParser(usage = "%prog [options] counts_file groups_file", description="LeafCutter differential splicing command line tool. Required inputs:\n <counts_file>: Intron usage counts file. Must be .txt or .txt.gz, output from clustering pipeline.\n <groups_file>: Two+K column file: 1. sample names (must match column names in counts_file), 2. groups (currently only two groups, i.e. pairwise, supported. Some samples in counts_file can be missing from this file, in which case they will not be included in the analysis. Additional columns can be used to specify confounders, e.g. batch/sex/age. Numeric columns will be treated as continuous, so use e.g. batch1, batch2, batch3 rather than 1, 2, 3 if you a categorical variable.",option_list=list(
  make_option(c("-o","--output_prefix"), default = "leafcutter_ds", help="The prefix for the two output files, <prefix>_cluster_significance.txt (containing test status, log likelihood ratio, degree of freedom, and p-value for each cluster) and <prefix>_effect_sizes.txt (containing the effect sizes for each intron)  [default %default]"),
  make_option(c("-s","--max_cluster_size"), default=Inf, help="Don't test clusters with more introns than this [default %default]"), 
  make_option(c("-i","--min_samples_per_intron"), default=5, help="Ignore introns used (i.e. at least one supporting read) in fewer than n samples [default %default]") , 
  make_option(c("-g","--min_samples_per_group"), default=3, help="Require this many samples in each group to have at least min_coverage reads [default %default]"), 
  make_option(c("-c","--min_coverage"), default=20, help="Require min_samples_per_group samples in each group to have at least this many reads [default %default]"), 
  make_option(c("-t","--timeout"), default=30, help="Maximum time (in seconds) allowed for a single optimization run [default %default]"),
  make_option(c("-p","--num_threads"), default=1, help="Number of threads to use [default %default]"),
  make_option(c("-e","--exon_file"), default=NULL, help="File defining known exons, example in data/gencode19_exons.txt.gz. Columns should be chr, start, end, strand, gene_name. Optional, only just to label the clusters."),
  make_option(c("--init"), default="smart", help="One of 'smart' (default) or 'random'."), 
  make_option(c("--seed"), default=12345, help="Random seed if using random initialization."))),
  positional_arguments = 2)

opt=arguments$opt
counts_file=arguments$args[1]
groups_file=arguments$args[2]

cat("Loading counts from",counts_file,"\n")
if (!file.exists(counts_file)) stop("File ",counts_file," does not exist")
counts=read.table(counts_file, header=T, check.names = F)

cat("Loading metadata from",groups_file,"\n")
if (!file.exists(groups_file)) stop("File ",groups_file," does not exist")
meta=read.table(groups_file, header=F, stringsAsFactors = F)
colnames(meta)[1:2]=c("sample","group")

counts=counts[,meta$sample]

group_names=unique(meta$group) # keep order from groups_file unless numeric
if (is.numeric(meta$group)) group_names=sort(group_names)
meta$group=factor(meta$group, group_names)

stopifnot(length(group_names)==2)

cat("Encoding as",group_names[1],"=0,",group_names[2],"=1\n")
numeric_x=as.numeric(meta$group)-1

confounders=NULL
if (ncol(meta)>2) {
    confounders=meta[,3:ncol(meta),drop=F]
    # scale continuous confounders
    for (i in seq_len(ncol(confounders)))
        if (is.numeric(confounders[,i]))
            confounders[,i]=scale(confounders[,i])
    # convert factors to one-of-K encoding
    confounders=model.matrix( ~., data=confounders )
    confounders=confounders[,2:ncol(confounders),drop=F] # remove intercept
}

minimum_group_size=min(sum(numeric_x==0),sum(numeric_x==1))
if (minimum_group_size < opt$min_samples_per_intron)
  stop("The number of samples in the smallest group is less than min_samples_per_intron, which means no clusters are testable. You can reduce min_samples_per_intron using the -i option, but note that we have only carefully checked the calibration of leafcutter p-values down to n=4 samples per group.")
if (minimum_group_size < opt$min_samples_per_group)
  stop("The number of samples in the smallest group is less than min_samples_per_group, which means no clusters are testable. You can reduce min_samples_per_intron using the -g option, but note that we have only carefully checked the calibration of leafcutter p-values down to n=4 samples per group.")

require(doMC)
registerDoMC(opt$num_threads)

cat("Settings:\n")
print(opt)

cat("Running differential splicing analysis...\n")
results <- differential_splicing(counts, numeric_x, confounders=confounders, max_cluster_size=opt$max_cluster_size, min_samples_per_intron=opt$min_samples_per_intron, min_samples_per_group=opt$min_samples_per_group, min_coverage=opt$min_coverage, timeout=opt$timeout, init=opt$init, seed=opt$seed ) 

cat("Saving results...\n")

# Make cluster table
cluster_table          = cluster_results_table(results)
cluster_table$cluster  = add_chr(cluster_table$cluster)

# Add gene names to clusters if an exon file is available
if (!is.null(opt$exon_file)) {
  cat("Loading exons from",opt$exon_file,"\n")
  if (file.exists(opt$exon_file)) {
     tryCatch( {
          exons_table     = read.table(opt$exon_file, header=T, stringsAsFactors = F)
          intron_meta     = get_intron_meta(rownames(counts))
          exons_table$chr = add_chr(exons_table$chr)
          intron_meta$chr = add_chr(intron_meta$chr)
          clu_gene_map    = map_clusters_to_genes(intron_meta, exons_table)
          cluster_table   = merge(cluster_table, clu_gene_map, by.x="cluster", by.y="clu", all.x=TRUE)
     }, error=function(err) warning(as.character(err)) ) # ignore errors here
  } else warning("File ",opt$exon_file," does not exist")
} else cat("No exon_file provided.\n")
write.table( cluster_table, paste0(opt$output_prefix,"_cluster_significance.txt"), quote=F, sep="\t", row.names = F)

# Write effect size table
effect_size_table                = leaf_cutter_effect_sizes(results)
colnames(effect_size_table)[3:4] = group_names
effect_size_table$intron = add_chr(effect_size_table$intron)
write.table( effect_size_table, paste0(opt$output_prefix,"_effect_sizes.txt"), quote=F, col.names = T, row.names = F, sep="\t")

# Save RData image
# save.image(paste0(opt$output_prefix,"_cluster_significance.RData"));
cat("All done, exiting\n")
