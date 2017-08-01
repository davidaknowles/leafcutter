#!/usr/bin/env Rscript
library(shiny, quietly = TRUE)
library(dplyr, quietly = TRUE)
library(ggplot2, quietly = TRUE)
library(DT, quietly = TRUE)
library(leafcutter,quietly = TRUE)
library(reshape2, quietly = TRUE)
library(gridExtra, quietly = TRUE)
library(intervals, quietly = TRUE) # needed for pretty strand arrow placement
library(ggplot2, quietly = TRUE)
library(foreach, quietly = TRUE)
library(shinycssloaders, quietly = TRUE)
library(optparse)
#data.table is required
exons_table <- NULL
resultsData <- NULL

# get Rdata file from option
opt <- parse_args(
  OptionParser(
    usage = "%prog results.Rdata",
    option_list=list()
  ), positional_arguments = c(0,1)
)

resultsData <- opt$args
if(length(resultsData) > 0){
  print(paste0("results file is ",resultsData) )
}else{
  resultsData <- NULL
}

if(!is.null(resultsData)){
  print("Loading in results")
  load(resultsData)

  # soon deprecated - exons_table will be loaded in with Rdata
  if(is.null(exons_table) ){
    if( species == "human" ){
      gencode_exons <- "data/gencode_hg38_all_exons.txt"
    }
    if( species == "mouse" ){
      gencode_exons <- "data/gencode_mm10_all_exons.txt"
    }
    print("reading in exons")
    exons_table <- as.data.frame(data.table::fread(gencode_exons))
  }
}else{
  print("no results found! load with default dataset")
  # remove variable from environment
  # server.R will then check if it exists
  rm(resultsData)
}

shiny::runApp( launch.browser=TRUE) 
