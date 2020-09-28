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
results_filename <- NULL

# get Rdata file from option
opt <- parse_args(
  OptionParser(
    usage = "%prog results.Rdata",
    option_list=list()
  ), positional_arguments = c(0,1)
)

results_filename <- opt$args
if(length(results_filename) == 0){
  print("No results found! Using default dataset example/Brain_vs_Heart_results.Rdata")
  results_filename <- "example/Brain_vs_Heart_results.Rdata"
}

print(paste0("Loading results from ",results_filename))
load(results_filename)

# option for Windows Subsystem for Linux
if (getOption("browser") == ""){
print("setting browser option to firefox")
options(browser = "firefox")
}

shiny::runApp( launch.browser=TRUE ) 
