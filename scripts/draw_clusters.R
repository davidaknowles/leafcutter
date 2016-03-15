#!/usr/bin/env Rscript
library(optparse)
library(leafcutter)

# e.g. Rscript ../scripts/draw_clusters.R brain_liver.txt.gz chrX:clu_8,chrX:clu_9 -g brain_liver_meta.txt -e ../leafcutter/data/gencode19_exons.txt.gz

arguments <- parse_args(OptionParser(usage = "%prog [options] counts_file cluster_name", description="LeafCutter cluster drawing tool. Required inputs:\n <counts_file>: Intron usage counts file. Must be .txt or .txt.gz, output from clustering pipeline.\n <cluster_name>: e.g. chr21:clu20. Multiple clusters may be separated by commas (no spaces!)",
option_list=list(
  make_option(c("-o","--output"), default = "cluster.pdf", help="The output file  [default %default]"),
  make_option(c("-e","--exon_file"), default=NULL, help="File defining known exons, example in data/gencode19_exons.txt.gz. Columns should be chr, start, end, strand, gene_name."),
  make_option(c("-g","--groups_file"), default=NULL, help="Two column file: 1. sample names (must match column names in counts_file), 2. groups. Some samples in counts_file can be missing from this file, in which case they will not be drawn."))), 
  positional_arguments = 2)

opt=arguments$opt
counts_file=arguments$args[1]
cluster_names=strsplit(arguments$args[2], ",")[[1]]

print(cluster_names)

cat("Loading counts from",counts_file,"\n")
if (!file.exists(counts_file)) stop("File ",counts_file," does not exist")
counts=read.table(counts_file, header=T)

meta=if (!is.null(opt$groups_file)) {
  cat("Loading metadata from",opt$groups_file,"\n")
  if (!file.exists(opt$groups_file)) stop("File ",opt$groups_file," does not exist")
  setNames( read.table(opt$groups_file, header=F, stringsAsFactors = F), c("sample","group") )
} else {
  cat("No groups provided.\n")
  NULL
}

exon_table=if (!is.null(opt$exon_file)) {
  cat("Loading exons from",opt$exon_file,"\n")
  if (!file.exists(opt$exon_file)) stop("File ",opt$exon_file," does not exist")
  read.table(opt$exon_file, header=T, stringsAsFactors = F)
} else {
  cat("No exon_file provided.\n")
  NULL
}

if (!is.null(meta)) {
  counts=counts[,meta$sample]
  meta$group=as.factor(meta$group)
  group_names=levels(meta$group)
}

introns=leafcutter:::get_intron_meta(rownames(counts))
cluster_ids=paste(introns$chr,introns$clu,sep = ":")

pdf(opt$output, width=8, height=8)
y=t(counts[ cluster_ids %in% cluster_names, ])
make_differential_splicing_plot(y, if (is.null(meta)) (numeric(nrow(y))+1) else meta$group, exons_table=exon_table)
dev.off()

cat("All done, exiting\n")
