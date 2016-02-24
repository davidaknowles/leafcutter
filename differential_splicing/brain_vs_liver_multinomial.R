# SETWD TO THIS FILE'S LOCATION FIRST

require(doMC)

source("../differential_splicing/differential_splicing.R")
source("../junction_plot/junction_plot.R")

registerDoMC(if (parallel::detectCores()==4) 7 else parallel::detectCores())

load("../example_data/brain_liver.RData") -> a

meta=as.data.frame(meta)
x=meta$tissue=="brain"

#res_file="../processed_data/brain_vs_liver_mult.RData"

#if (!file.exists(res_file)) {
  results=differential_splicing(numers, x, max_cluster_size = 5)
#  save(results, file=res_file)
#} else load(res_file)

cluster_table=cluster_results_table(results)
pqplot(cluster_table$p)

intron_effect_sizes=leaf_cutter_effect_sizes(results)

top_clus=rownames(cluster_table)[ order(cluster_table$p)[1:10] ]
counter=1
toWrite=list()
for (clu in top_clus) {  
  #fn=paste0("../brain_vs_heart_examples/brain_vs_heart_",counter,".png")
 # png(fn,width = 7,height=7,units = "in", res = 150)
  y=t(numers[ cluster_ids==clu, ])
  exons=make_differential_splicing_plot(y, meta$tissue)
  counter = counter+1
  #dev.off())
}