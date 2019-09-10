#!/usr/bin/env Rscript
#options(echo=TRUE)
#### Jack Humphrey 2017
### annotate the output of differential splicing
## and prepare for visualisation

library(optparse)

option_parser=OptionParser(
  usage="%prog [options] <name>_perind_numers.counts.gz <name>_cluster_significance.txt <name>_effect_sizes.txt annotation_code \nThe annotation_code should be something like annotation_codes/gencode_hg19/gencode_hg19.",
  option_list=list(
    make_option( c("-o","--output"), default="leafviz.RData", help="The output file that will be created ready for loading by run_leafviz.R [%default]"),
    make_option( c("-m","--meta_data_file"), default=NULL, help="The support file used in the differential splicing analysis. Columns should be file name and condition"),
  make_option( c("-f","--FDR"), default=0.05, help = "the adjusted p value threshold to use [%default]"),
    make_option( c("-c","--code"), default="leafcutter_ds", help = "A name for this analysis (will be available in leafviz through the Summary tab). [%default]"))
)

parsed_args <- parse_args(option_parser,  positional_arguments = 4)

counts_file <- parsed_args$args[1]
cluster_significance_file <- parsed_args$args[2]
effect.sizes.file <- parsed_args$args[3]
annotation_code <- parsed_args$args[4]

code <- parsed_args$options$code
results_file = parsed_args$options$output
groups_file <- parsed_args$options$meta_data_file
FDR_limit <- parsed_args$options$FDR

require(leafcutter)
library(data.table)
library(stringr)
library(dplyr)
library(magrittr)

cat("Preparing for visualisation\n")

# annotation
exon_file <- paste0(annotation_code, "_all_exons.txt.gz")
all_introns <- paste0(annotation_code,"_all_introns.bed.gz" )
threeprime_file <- paste0( annotation_code,"_threeprime.bed.gz")
fiveprime_file <- paste0( annotation_code,"_fiveprime.bed.gz")

# CHECK DEPENDENCY FILES
pass <- TRUE
errorMessage <- c()
for( file in
  c(effect.sizes.file, cluster_significance_file, counts_file, groups_file,
    all_introns, threeprime_file, fiveprime_file, exon_file
    )){
  if( !file.exists(file)){
    pass <- FALSE
    errorMessage <- c(errorMessage,  paste0(file, " does not exist\n") )
  }
}
if(!pass){
  stop(errorMessage)
}
cat("Results to be saved in:",results_file, "\n")
cat("Using annotation at:", annotation_code,"\n")

cat("Loading counts from",counts_file,"\n")
counts <- read.table(counts_file, check.names=FALSE)

if(file.exists(groups_file)){ # can we run without this?
  cat("Loading metadata from",groups_file,"\n")
  meta <- read.table(groups_file, header=FALSE, stringsAsFactors = FALSE)
  colnames(meta)=c("sample","group")
  # name covariates
  if( ncol(meta) > 2){
    colnames(meta)[3:ncol(meta)] <- paste0("covariate", 1:(ncol(meta) - 2)  )
  }
  sample_table <- data.frame( group = names(table(meta$group) ), count = as.vector(table(meta$group)) )
}

# exon table no longer used for anything - just saved with the Rdata object at the end
exons_table=if (!is.null( exon_file )) {
  cat("Loading exons from",exon_file,"\n")
  #read_table(exon_file)
  as.data.frame(fread(paste("zless",exon_file)), data.table=F )
} else {
  cat("No exon_file provided.\n")
  NULL
}

# in case the exons table lacks "chr" on the chr column
exons_table$chr <- leafcutter::add_chr(exons_table$chr)

effectSizes <- fread(effect.sizes.file, data.table=F )
effectSizesSplit <-  as.data.frame(str_split_fixed(effectSizes$intron, ":", 4), stringsAsFactors = FALSE )
names(effectSizesSplit) <- c("chr","start","end","clusterID")

effectSizes <- cbind( effectSizes, effectSizesSplit)
effectSizes$cluster <- paste(effectSizesSplit$chr, effectSizesSplit$clusterID, sep = ":")

results <- fread(cluster_significance_file, data.table=F )
results$FDR <- p.adjust( results$p, method = "fdr")

# If there were no significant results, stop here and provide feedback to the user by printing
# an error message indicating that no significant clusters were found at this FDR threshold
if( !any(results$FDR < FDR_limit, na.rm=T) ){
   stop("No significant clusters found\n");
}

# Gather introns meeting the FDR threshold
all.introns <- merge(x = results, y = effectSizes, by = "cluster")

if( nrow(all.introns) == 0 ){
  stop("Merging the per-cluster results with the per-junction effect sizes produces an empty table. Please check your input files.")
}

all.introns <- all.introns[ order(all.introns$FDR),]

all.introns <- subset( all.introns, FDR <= FDR_limit )
all.introns$start <- as.numeric(all.introns$start)
all.introns$end <- as.numeric(all.introns$end)

# for introns, 5' splice sites and 3 splice sites:
# add "chr" to chrom name if needed
## intersect with list of junctions
all.junctions <- dplyr::select(all.introns, chr, start, end, clusterID)

intron_db <- fread(paste0("zcat < ", all_introns), data.table = FALSE)
colnames(intron_db)[1:4]=c("chr","start","end","gene")
intron_db$chr <- leafcutter::add_chr(intron_db$chr)
all.introns_intersect = all.junctions %>%
  left_join(intron_db, by=c("chr","start","end"))

threeprime_db <- fread(paste0("zcat < ", threeprime_file), data.table = FALSE)
colnames(threeprime_db)[1:7]=c("chr","start","end","gene","gene_id","strand","transcript")
threeprime_db$chr <- leafcutter::add_chr(intron_db$chr)
threeprime_intersect = all.junctions %>%
  select(chr, clusterID, start=end) %>%
  left_join(threeprime_db, by=c("chr","start"))

fiveprime_db <- fread(paste0("zcat < ", fiveprime_file), data.table = FALSE)
colnames(fiveprime_db)[1:7]=c("chr","start","end","gene","gene_id","strand","transcript")
fiveprime_db$chr <- leafcutter::add_chr(fiveprime_db$chr)
fiveprime_intersect =  all.junctions %>%
  select(chr, clusterID, start) %>%
  left_join(fiveprime_db, by=c("chr","start"))

# now I have two lists of splice site annotation
# for testing
print("Annotating junctions")

verdict.list <- list()
coord.list <- list()
gene.list <- list()
ensemblID.list <- list()
transcripts.list <- list()
constitutive.list <- list()
classification.list <- list()

clusters <- unique( all.introns$clusterID )
for( clu in clusters ){

  # for each intron in the cluster, check for coverage of both
  # output a vector of string descriptions
  cluster <- all.introns %>% filter( clusterID == clu )

  # first subset the intersected files to speed up later query - this uses the data.tables method
  fprimeClu <- fiveprime_intersect %>% filter( clusterID == clu )
  tprimeClu <- threeprime_intersect %>% filter( clusterID == clu )
  bothSSClu <- all.introns_intersect %>% filter( clusterID == clu )

  # for each intron in the cluster:
  #   create vector of overlapping splice sites, indexed by the row of the intersect
  # five prime splice sites
  fprime=cluster %>% left_join(fprimeClu, by=c("chr","start"))

  # three prime splice sites
  tprime=cluster %>% left_join(tprimeClu, by=c("chr"="chr","end"="start"))

  # both splice sites
  bothSS=cluster %>% left_join(bothSSClu, by=c("chr","start","end"))

  # find gene and ensemblID by the most represented gene among all the splice sites - lazy
  cluster_gene <- names(sort(table(c(tprime$gene,fprime$gene)), decreasing = TRUE ))[1]

    # if no cluster gene found then leave as "."
  if( is.null(cluster_gene) ){
    cluster_gene <- "."
  }

  gene_strand <- NA
  if( cluster_gene != "." ){
    # get strand the same way - would prefer to use the strand of the junction
    strands <- c(tprime$strand, fprime$strand)
    # hope that all junctions align to the same gene on the same strand
    gene_strand <- unique( strands[ strands != "." & !is.na(strands) ])
    if( all(is.na(gene_strand)) | length(gene_strand) != 1 ){
      gene_strand <- NA
    }
  }

  # do the same for EnsemblID
  cluster_ensemblIDs <- names(sort(table( c(tprime$gene_id,fprime$gene_id)), decreasing = TRUE ))
  cluster_ensemblID <- cluster_ensemblIDs[ cluster_ensemblIDs != "." ][1]
  if( length( cluster_ensemblID ) == 0 ){
    cluster_ensemblID == "."
  }

  verdict <- c()
  coord <- c()
  gene <- c()
  ensemblID <- c()
  transcripts <- list()

  for( intron in 1:nrow(cluster) ){
    coord[intron] <- paste(cluster[intron,]$chr,cluster[intron,]$start, cluster[intron,]$end )

    gene[intron] <- cluster_gene
    ensemblID[intron] <- cluster_ensemblID

    fprime_intron=cluster[intron,] %>% left_join(fprime, by=c("chr","start"))
    tprime_intron=cluster[intron,] %>% left_join(tprime, by=c("chr","end"))
    bothSS_intron=cluster[intron,] %>% left_join(bothSSClu, by=c("chr","start","end"))

    # for each intron create vector of all transcripts that contain both splice sites
    transcripts[[intron]] <- intersect( tprime_intron$transcript,fprime_intron$transcript )

    verdict[intron] <- "error"

    unknown_3p=all( is.na(tprime_intron$gene) )
    unknown_5p=all( is.na(fprime_intron$gene) )

    if (is.na(gene_strand)) {
      verdict[intron] <- "unknown_strand"
    } else {
      if( all( is.na(tprime_intron$gene )) & all( is.na(fprime_intron$gene ) ) ){
        verdict[intron] <- "cryptic_unanchored"
      }
      if( (all( is.na(tprime_intron$gene )) & all( !is.na(fprime_intron$gene ) ) & all(gene_strand == "+") ) |
        ( all( is.na(fprime_intron$gene )) & all( !is.na(tprime_intron$gene ) ) & all(gene_strand == "-") )
      ){ verdict[intron] <- "cryptic_threeprime"
      }
      if(
        ( all( !is.na(tprime_intron$gene )) & all( is.na(fprime_intron$gene ) ) & all(gene_strand == "+") ) |
        ( all( !is.na(fprime_intron$gene )) & all( is.na(tprime_intron$gene ) ) & all(gene_strand == "-") )
      ){ verdict[intron] <- "cryptic_fiveprime"
      }
      if( is.na(gene_strand) & ( all( !is.na(tprime_intron$gene )) | all( !is.na(fprime_intron$gene ) ) ) ){
        verdict[intron] <- "cryptic"
      }
      if( # if both splice sites are annotated
        all( !is.na(tprime_intron$gene ) ) & all( !is.na(fprime_intron$gene ) )
      ){
        # test if the splice sites are paired in a known intron
        if( all( !is.na(bothSS_intron$gene )) ){
          verdict[intron] <- "annotated"
        }else{ # both are annotated but never in the same junction
          verdict[intron] <- "novel annotated pair"
        }
      }
    }
    verdict.list[[clu]] <- verdict
    coord.list[[clu]] <- coord
    gene.list[[clu]] <- gene
    ensemblID.list[[clu]] <- ensemblID
    #transcripts.list[[clu]] <- transcripts

    # once all the transcripts for all the introns are found, go back and work out how many constitutive each junction is. Does the junction appear in every transcript?

    if( intron == nrow(cluster)){ # only on final intron
      all_transcripts <- unique( unlist( transcripts ) )
      # remove "." - non-existent transcripts
      all_transcripts <- all_transcripts[ all_transcripts != "." ]

      constitutive <- lapply( transcripts, FUN = function(x) {
        # for each intron how many transcripts is it seen in?
        x <- x[ x != "." ]
        length(x) / length( all_transcripts)

        })

      constitutive.list[[clu]] <- constitutive

      # collapse all.introns transcripts for each intron into a single string
      transcripts.list[[clu]] <- lapply(transcripts, FUN = function(x) paste( x, collapse = "+" ) )

    }

  }

  # predicting the event type from the shape of the junctions
  #print(clu)

  if( nrow(cluster) != 3){
    classification.list[[clu]] <- "."
    next
  }else{
    classification.list[[clu]] <- "."

    tab <- select(cluster, start, end)

    # the junctions are sorted by start and end coordinates

    # check for the presence of a junction that spans the entire length of the cluster
    if( !any(  which( tab$start == min(tab$start) ) %in% which( tab$end == max(tab$end) )  ) ){
      classification.list[[clu]] <- "."
      next
    }

    # therefore for a cassette exon arrangement the longest junction always comes second
    if( which( tab$start ==  min(tab$start) & tab$end == max(tab$end ) ) != 2 ){
     classification.list[[clu]] <- "."
     next
    }

    # now we know that junction 2 is the parent, junction 1 is the left most child and junction 3 is the right most
    # check that the end of junction 1 comes before the start of junction 3

    if( tab[1,"end"] > tab[3,"start"] ){
      classification.list[[clu]] <- "."
      next
    }

    # double check the starts and ends
    if( tab[1, "start"] != tab[2,"start"] | tab[3,"end"] != tab[2,"end"] ){
      classification.list[[clu]] <- "."
      next
    }

    # work out direction of change
    if( cluster[1, "deltapsi"] > 0 & cluster[3, "deltapsi"] > 0 & cluster[2,"deltapsi"] < 0){
      classification.list[[clu]] <- "cassette exon - increased"
    }
    if( cluster[1, "deltapsi"] < 0 & cluster[3, "deltapsi"] < 0 & cluster[2,"deltapsi"] > 0){
      classification.list[[clu]] <- "cassette exon - decreased"
    }

    # work out annotation status
    if( all( verdict.list[[clu]] == "annotated") ){
      classification.list[[clu]] <- paste0( classification.list[[clu]], " - annotated")
    }

    if( verdict.list[[clu]][2] == "annotated" & verdict.list[[clu]][1] != "annotated" & verdict.list[[clu]][3] != "annotated"  ){
      classification.list[[clu]] <- paste0( classification.list[[clu]], " - cryptic")
    }

    if( verdict.list[[clu]][2] == "novel annotated pair" & verdict.list[[clu]][1] == "annotated" & verdict.list[[clu]][3] == "annotated"  ){
      classification.list[[clu]] <- paste0( classification.list[[clu]], " - skiptic")
    }

  }

}

#save.image(file = "test.Rdata")

print("Preparing results")

# match the lists together
all.introns$verdict <- unlist(verdict.list)[ match( paste( all.introns$chr, all.introns$start, all.introns$end ), unlist(coord.list)) ]

all.introns$gene <- unlist(gene.list)[ match( paste( all.introns$chr, all.introns$start, all.introns$end ), unlist(coord.list)) ]

all.introns$ensemblID <- unlist(ensemblID.list)[ match( paste( all.introns$chr, all.introns$start, all.introns$end ), unlist(coord.list)) ]

all.introns$transcripts <- unlist( transcripts.list )[ match( paste( all.introns$chr, all.introns$start, all.introns$end ), unlist(coord.list)) ]

all.introns$constitutive.score <-  unlist( constitutive.list )[ match( paste( all.introns$chr, all.introns$start, all.introns$end ), unlist(coord.list)) ]

#all.introns$prediction <-  unlist( classification.list )[ match( paste( all.introns$chr, all.introns$start, all.introns$end ), unlist(coord.list)) ]

# replace NA values/missing transcripts with "."
all.introns %<>% mutate( gene=ifelse(is.na(gene), ".", gene),
                         ensemblID=ifelse(is.na(ensemblID), ".", ensemblID),
                         transcripts=ifelse(transcripts == "", ".", transcripts),
                         constitutive.score=signif(constitutive.score, digits = 2))

# prepare results
results$clusterID <- str_split_fixed(results$cluster, ":", 2)[,2]
results$N <- results$df + 1
sig <- subset(results, FDR < FDR_limit)
sig$clusterID <- str_split_fixed(sig$cluster, ":", 2)[,2]

all.clusters <- lapply(sig$clusterID, FUN = function(clu){
  cluster <- all.introns[ all.introns$clusterID == clu, ]
  chr <- unique( cluster$chr )[1] # this should always be one number
  start <- min( cluster$start )
  end <- max( cluster$end )
  # get most common gene name that is not "."
  gene <- names( sort( table( unique(cluster$gene) ), decreasing = TRUE ) )[1]
  ensemblID <- names( sort( table( unique(cluster$ensemblID) ), decreasing = TRUE ) )[1]
  annotation <- "annotated"
  if( any(grepl( "cryptic", cluster$verdict)) | any( grepl("novel annotated pair", cluster$verdict)) ){
    annotation <- "cryptic"
  }
  return(
    data.frame(
      clusterID = clu,
      chr = chr,
      start = start,
      end = end,
      gene = gene,
      ensemblID = ensemblID,
      annotation = annotation ) )
  })
all.clusters <- do.call( what = rbind, args = all.clusters)

all.clusters$FDR  <- results$FDR[ match( all.clusters$clusterID, results$clusterID)]
all.clusters$FDR <- signif( all.clusters$FDR, digits = 3)
all.clusters$N  <- results$N[ match( all.clusters$clusterID, results$clusterID)]

# add classification
all.clusters$verdict <- unlist(classification.list)[ match(all.clusters$clusterID, names(classification.list))]

# prepare for PCA
counts <- counts[,meta$sample]
print( "converting counts to ratios")
# create per cluster ratios from counts
ratios <- counts %>%
  mutate(clu = str_split_fixed(rownames(counts), ":", 4)[,4]) %>%
  group_by(clu) %>%
  mutate_all( funs( ./sum(.) ) ) %>%
  ungroup() %>%
  as.data.frame() %>%
  magrittr::set_rownames(rownames(counts)) %>%
  select(-clu)
ratios <- ratios[rowMeans(is.na(ratios)) <= 0.4,,drop=FALSE ]
row_means <- rowMeans(ratios, na.rm = TRUE)
row_means_outer <- outer(row_means, rep(1,ncol(ratios)))
ratios[is.na(ratios)] <- row_means_outer[is.na(ratios)]


meta$group <- as.factor(meta$group)

make_pca <- function(counts,meta){
  dev <- apply( counts, MAR = 1, FUN = sd )
  # remove rows with 0 variance
  counts <- counts[ dev != 0, ]
  pca <- prcomp( t(counts), scale = TRUE )
  importance <- signif( summary(pca)$importance[2,], digits = 2) * 100
  pca <- as.data.frame(pca$x)
  pca$sample <- row.names(pca)
  pca <- merge(pca,meta, by = "sample")
  row.names(pca) <- pca$sample
  pca$sample <- NULL

  return(list( pca, importance) )
}
print("creating PCA")
pca <- make_pca(ratios,meta)

# sort out clusters table
fix_clusters <- function(clusters){
  clusters$FDR <- signif( clusters$FDR, digits = 3)
  clusters$coord <- paste0( clusters$chr, ":", clusters$start, "-", clusters$end)
  clusters <- clusters[ order(clusters$FDR, decreasing = FALSE),]
  # removed ensemblID - this could be an option?
  #clusters <- select( clusters, clusterID, N, coord, gene, annotation, FDR, verdict)
  clusters <- select(
    clusters,
    clusterID,
    N,
    coord,
    gene,
    annotation,
    FDR)
  clusters$gene <- paste0("<i>",clusters$gene,"</i>")
return(clusters)
}

# use on all.introns
fix_introns <- function(introns){
  introns <- select(introns,
    clusterID,
    gene,
    ensemblID,
    chr,
    start,
    end,
    verdict,
    deltapsi,
    #constitutive.score,
    transcripts)
  #introns$constitutive.score <- signif(introns$constitutive.score, digits = 3)
  introns$deltapsi<- round(introns$deltapsi, digits = 3)
return(introns)
}


cluster_summary <- function(clusters){
  summary <- data.frame(
              Results = c(
                paste0("Number of differentially spliced clusters at FDR = ", FDR_limit) ,
                        "Fully annotated",
                        "Contain unannotated junctions"),
                n = c( nrow(clusters),
                       nrow( clusters[ clusters$annotation == "annotated", ]),
                       nrow( clusters[ clusters$annotation == "cryptic", ])
                       )
                )
  return(summary)
}

intron_summary <- function(all.introns){
    summary <- data.frame(
                Results = c( "Number of fully annotated junctions",
                             "Number of junctions with cryptic 5' splice site",
                              "Number of junctions with cryptic 3' splice site",
                              "Number of junctions with two cryptic splice sites",
                              "Number of novel junctions that connect two annotated splice sites"),

                  n = c( nrow(all.introns[ all.introns$verdict == "annotated",]),
                         nrow(all.introns[ all.introns$verdict == "cryptic_fiveprime",]),
                         nrow(all.introns[ all.introns$verdict == "cryptic_threeprime",]),
                         nrow(all.introns[ all.introns$verdict == "cryptic_unanchored",]),
                         nrow(all.introns[ all.introns$verdict == "novel annotated pair",])
                  )
                )
    return( summary )
}

# create all the objects for visualisation
clusters <- fix_clusters(all.clusters)
introns <- fix_introns(all.introns)
intron_summary <- intron_summary(all.introns)
cluster_summary <- cluster_summary(all.clusters)
introns_to_plot <- get_intron_meta(rownames(counts))
cluster_ids <- introns_to_plot$clu

# save all the objects needed by Leafcutter viz into single Rdata file
# include the mode variable
save( introns,
      clusters,
      counts,
      meta,
      exons_table,
      pca,
      intron_summary,
      cluster_summary,
      introns_to_plot,
      cluster_ids,
      sample_table,
      annotation_code,
      code,
      file = results_file
)

