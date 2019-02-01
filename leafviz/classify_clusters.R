
# leafcutter cluster prediction
# Jack Humphrey
# 2017 - 2019

# input: a LeafViz RData file

library(optparse)
library(dplyr, quietly = TRUE,warn.conflicts = FALSE)

option_parser=OptionParser(
  usage="%prog [options] <your analysis>_leafviz.Rdata",
  option_list=list(
    make_option( c("-o","--output"), 
                 default="clusters", 
                 help="The name and path of where to output the classified clusters [%default]")#,
   # make_option( c("-m","--meta_data_file"), default=NULL, help="The support file used in the differential splicing analysis. Columns should be file name and condition"),
  #make_option( c("-f","--FDR"), default=0.05, help = "the adjusted p value threshold to use [%default]"),
  #  make_option( c("-c","--code"), default="leafcutter_ds", help = "A name for this analysis (will be available in leafviz through the Summary tab). [%default]"))
  )
)

parsed_args <- parse_args(option_parser,  positional_arguments = 1)

leafviz_data <- parsed_args$args[1]

#code <- parsed_args$options$code
results_base = parsed_args$options$output
#groups_file <- parsed_args$options$meta_data_file
#FDR_limit <- parsed_args$options$FDR

#require(leafcutter)
#library(data.table)

#library(stringr)
#library(dplyr)
#library(magrittr)

# load in data

message(paste0("loading data from ", leafviz_data))

load(leafviz_data)


classifyClusters <- function(clu, full_output = FALSE){
  
  # retrive cluster 
  cluster <- filter(introns, clusterID == clu )

  ## ----------------------------------------------------
  ## CASSETTE EXON PREDICTION AND CLASSIFICATION --------
  
  # most clusters will not conform so output a null result
  # usually unnecessary to include this though
  if( full_output == TRUE ){ 
    null_class <- 
      data.frame(
      clusterID = clu,
      verdict = "other",
      direction = ".",
      junctions_annotated = ".",
      exon_annotated = "."
    )
  }else{
    null_class <- NULL
  }

  # only works on 3 junction clusters 
  if( nrow(cluster) != 3){ 
    return(null_class)
  }
  
  # the junctions are sorted by start and end coordinates
  tab <- select(cluster, start, end)
  

  ## TOPOLOGY TESTS -------------
  
  # check for the presence of a junction that spans the entire length of the cluster
  if( !any(  which( tab$start == min(tab$start) ) %in% which( tab$end == max(tab$end) )  ) ){
    return(null_class)
    }
  
  # therefore for a cassette exon arrangement the longest junction always comes second 
  if( which( tab$start ==  min(tab$start) & tab$end == max(tab$end ) ) != 2 ){
    return(null_class)
    }
  
  # now we know that junction 2 is the parent, junction 1 is the left most child and junction 3 is the right most
  # check that the end of junction 1 comes before the start of junction 3
  
  if( tab[1,"end"] > tab[3,"start"] ){
    return(null_class)
  }
  
  # double check the starts and ends
  if( tab[1, "start"] != tab[2,"start"] | tab[3,"end"] != tab[2,"end"] ){
    return(null_class)
  }
  
  # if these tests are passed then the variant is indeed a cassette exon
  class <- "cassette"
  
  # work out direction of change
  # some exons don't get classified because the deltapsi parameter doesn't match the direction of change.
  # but if the parent junction goes down and one of the two children goes up then that's enough for me.
  direction <- "unclear"
  
  if( cluster[2, "deltapsi"] < 0 & ( cluster[3, "deltapsi"] > 0 | cluster[1,"deltapsi"] > 0 )  ){
    direction <- "included"
  }
  
  if( cluster[2, "deltapsi"] > 0 & ( cluster[3, "deltapsi"] < 0 | cluster[1,"deltapsi"] < 0 ) ){
    direction <- "skipped"
  }
  
  
  # 
  # junction annotation 
  # provide separate columns for skipping and inclusion junctions
  # for inclusion junctions - if both are annotated then its annotated; 
  # if one is novel then its novel
  # if both are novel then its novel
  # for skipping junction - either novel or annotated
  anno_skip <- "novel"
  anno_include <- "novel"
  
  if( cluster$verdict[2] == "annotated" ){
  anno_skip <- "annotated"
  }else{
    anno_skip <- "novel"
  }
  if( cluster$verdict[1] == "annotated" & cluster$verdict[3] == "annotated" ){
    anno_include <- "annotated"
  }else{
    anno_include <- "novel"
  }
  
  
  # Exon annotation
  # if cassette is annotated, check that there is annotation to support an exon joining the two junctions
  anno_exon <- "novel"
  
  exon.start <- cluster$end[1]
  exon.end <- cluster$start[3]
  
  exon_overlap <- dplyr::filter(exons_table, start == exon.start, end == exon.end)
  
  if(nrow(exon_overlap) > 0){
    anno_exon <- "annotated"
  }
  
  


  # return verdict - what people want
  # if both junctions and exon are annotated it is a cassette exon with a direction
  # 
  verdict <- "complex"
  if( anno_skip == "annotated" & anno_include == "annotated" & anno_exon == "annotated" & direction != "unclear"){
    verdict = paste0("cassette:",direction)
  }
  if( anno_skip == "novel" & anno_include == "annotated" & anno_exon == "annotated" & direction != "unclear"){
    verdict = paste0("novelskip:",direction)
  }
  if( anno_skip == "annotated" & anno_include == "novel" & anno_exon == "novel" & direction != "unclear"){
    verdict = paste0("novelexon:",direction)
  }
  
  # paste together predictions
  # add gene and coords
  result <- 
    data.frame(
      clusterID = clu,
      gene = unique(cluster$gene),
      direction = direction,
      anno_skip = anno_skip,
      anno_include = anno_include,
      anno_exon = anno_exon,
      verdict = verdict
    )
  
  return(result)
}


message( "classifying clusters")

# remove purrr and tibble as dependencies
results <- 
  do.call(
    lapply(X = clusters$clusterID, 
           FUN = function(x){ 
             classifyClusters(x, full_output = FALSE) } ),
    what = "rbind"
  )


results_summary <-
  results %>%
  group_by(verdict) %>%
  summarise(n = n() ) %>%
  mutate( prop = n / sum(.$n)) %>%
  arrange(desc(n))

print(results_summary)

message("saving results")

write.table( results, 
             file = paste0(results_base,"_classifications.tsv"), 
             quote = FALSE, 
             sep = "\t",
             row.names = FALSE )

write.table( results_summary, 
             file = paste0(results_base,"_summary.tsv"), 
             quote = FALSE, 
             sep = "\t",
             row.names = FALSE )




