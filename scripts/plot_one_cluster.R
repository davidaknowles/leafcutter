#!/usr/bin/env Rscript
library(optparse)
library(leafcutter)

arguments <- parse_args(OptionParser(usage = "%prog [options] counts_file groups_file cluster_name", description="LeafCutter differential splicing plotting tool. Required inputs:\n <counts_file>: Intron usage counts file. Must be .txt or .txt.gz, output from clustering pipeline.\n <groups_file>: Two column file: 1. sample names (must match column names in counts_file), 2. groups (currently only two groups, i.e. pairwise, supported. Some samples in counts_file can be missing from this file, in which case they will not be included in the analysis.",option_list=list(
  make_option(c("-o","--output"), default = "leafcutter_plot.pdf", help="The output file  [default %default]"),
  make_option(c("-e","--exon_file"), default=NULL, help="File defining known exons, example in data/gencode19_exons.txt.gz. Columns should be chr, start, end, strand, gene_name."))), 
  positional_arguments = 3)

opt=arguments$opt
counts_file=arguments$args[1]
groups_file=arguments$args[2]
cluster_to_plot=arguments$args[3]

cat("Loading counts from",counts_file,"\n")
if (!file.exists(counts_file)) stop("File ",counts_file," does not exist")
counts=read.table(counts_file, header=T)

cat("Loading metadata from",groups_file,"\n")
if (!file.exists(groups_file)) stop("File ",groups_file," does not exist")
meta=read.table(groups_file, header=F, stringsAsFactors = F)
colnames(meta)=c("sample","group")

exon_table=if (!is.null(opt$exon_file)) {
  cat("Loading exons from",opt$exon_file,"\n")
  if (!file.exists(opt$exon_file)) stop("File ",opt$exon_file," does not exist")
  read.table(opt$exon_file, header=T, stringsAsFactors = F)
} else {
  cat("No exon_file provided.\n")
  NULL
}

counts=counts[,meta$sample]

meta$group=as.factor(meta$group)
group_names=levels(meta$group)

cat("Saving",cluster_to_plot,"plots to",opt$output,"\n")
introns=leafcutter:::get_intron_meta(rownames(counts))
cluster_ids=paste(introns$chr,introns$clu,sep = ":")

pdf(opt$output, width=8, height=8)
y=t(counts[ cluster_ids==cluster_to_plot, ])
make_differential_splicing_plot(y, meta$group, exons_table=exon_table)
dev.off()

cat("All done, exiting\n")
