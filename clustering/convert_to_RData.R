#!/usr/bin/env Rscript
args=commandArgs(trailingOnly = T)

numers=read.table(args[1], header=T, row.names=1, stringsAsFactors=F)

#ss=strsplit(colnames(numers),split="_")
#meta=do.call(rbind,ss)
#colnames(meta)=c("sex","tissue","id")

save(numers, file=args[2] ) 
