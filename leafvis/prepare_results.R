#!/usr/bin/env Rscript

#### Jack Humphrey 2017
### annotate the output of differential splicing
## and prepare for visualisation
library(optparse)

options(echo=TRUE)

opt <- parse_args(
  OptionParser(
    option_list=list(
      make_option( c("-i", "--iFolder"), help = "The folder where the LeafCutter differential splicing was run"),
      make_option( c("-o","--outFolder"), help = "The folder where the results object will be saved" ),
      make_option( "--support", default=NULL, help="The support file used in the differential splicing analysis. Columns should be file name and condition"),
      make_option( "--annotation_code", default=NULL, help = "a path to the annotation files of exons, introns and splice sites"),
      make_option( c("-f","--FDR"), default=0.05, help = "the adjusted p value threshold to use"),
      make_option( "--code", default=NULL, help = "the same dataset-specific code used throughout the pipeline"))
	 )
  )

outFolder <- opt$outFolder
iFolder <- opt$iFolder
species <- opt$species
groups_file <- opt$support
# counts_file <- opt$counts_file
annotation_code <- opt$annotation_code
code <- opt$code
FDR_limit <- opt$FDR




cat("Preparing for visualisation\n")

# results
effect.sizes.file <- paste0(iFolder,"/",code,"_ds_effect_sizes.txt") 
results.file <- paste0(iFolder, "/",code,"_ds_cluster_significance.txt")
counts_file <- paste0(iFolder, "/",code, "_perind_numers.counts.gz")
# annotation
exon_file <- paste0(annotation_code, "_all_exons.txt.gz")
all_introns <- paste0(annotation_code,"_all_introns.bed.gz" )
threeprime_file <- paste0( annotation_code,"_threeprime.bed.gz")
fiveprime_file <- paste0( annotation_code,"_fiveprime.bed.gz")


# CHECK DEPENDENCY FILES
pass <- TRUE
errorMessage <- c()

# Dependent on
# effects sizes and results file, not counts or groups

for( file in 
  c(effect.sizes.file, results.file, counts_file, groups_file, 
    all_introns, threeprime_file, fiveprime_file, exon_file
    )){
  if( !file.exists(file)){
    pass <- FALSE
    errorMessage <- c(errorMessage,  paste0(file, " does not exist\n") )
  }
}
if(!pass){
  #print(errorMessage)
  stop(errorMessage)
}
cat(paste("results to be saved in:",outFolder, "\n") )
cat(paste("using annotation found in:", annotation_code,"\n") )


library(data.table)
#library(leafcutter)
library(stringr)
library(dplyr)


# CHECK BEDTOOLS IS INSTALLED
if( !file.exists(system("which bedtools", intern=TRUE) )){
  stop("bedtools is not in your $PATH - have you installed it?")
}

if(file.exists(counts_file)){
  cat("Loading counts from",counts_file,"\n")
  counts <- read.table(counts_file)
}

if(file.exists(groups_file)){
  cat("Loading metadata from",groups_file,"\n")
  meta <- read.table(groups_file, header=F, stringsAsFactors = F)
  colnames(meta)=c("sample","group")

  sample_table <- data.frame( group = names(table(meta$group) ), count = as.vector(table(meta$group)) )
}
# exon table no longer used for anything - just saved with the Rdata object at the end

exons_table=if (!is.null( exon_file )) {
  cat("Loading exons from",exon_file,"\n")
  #read_table(exon_file)
  as.data.frame(fread(paste("zless",exon_file)) )
} else {
  cat("No exon_file provided.\n")
  NULL
}

resultsFolder <- outFolder
if( !dir.exists(resultsFolder)){
  dir.create(resultsFolder)
}


effectSizes <- fread(effect.sizes.file, stringsAsFactors = FALSE )
effectSizesSplit <-  as.data.frame(str_split_fixed(effectSizes$intron, ":", 4), stringsAsFactors = FALSE )
names(effectSizesSplit) <- c("chr","start","end","clusterID")

print(head(effectSizes))
print(head(effectSizesSplit))

effectSizes <- cbind( effectSizes, effectSizesSplit)
effectSizes$cluster <- paste(effectSizesSplit$chr, effectSizesSplit$clusterID, sep = ":")

results <- fread(results.file, stringsAsFactors = F)

results$FDR <- p.adjust( results$p, method = "fdr")

all.introns <- merge(x = results, y = effectSizes, by = "cluster")
all.introns <- all.introns[ order(all.introns$FDR),]

all.introns <- subset( all.introns, FDR <= FDR_limit )
all.introns$start <- as.numeric(all.introns$start)
all.introns$end <- as.numeric(all.introns$end)


# for each splice site write out a bed file  
all.junctions <- dplyr::select(all.introns, chr, start, end, clusterID)

all.fiveprime <- data.frame( chr = all.introns$chr,
                             start = all.introns$start,
                             end = as.numeric( as.character(all.introns$start) ) + 1,
                             clusterID = all.introns$clusterID)
all.threeprime <- data.frame( chr = all.introns$chr,
                             start = all.introns$end,
                             end = as.numeric( as.character(all.introns$end) ) + 1,
                             clusterID = all.introns$clusterID)
all.file <- paste0(resultsFolder, "/all_junctions.bed")
all.fiveprime.file <- paste0(resultsFolder, "/all_fiveprime.bed")
all.threeprime.file <- paste0(resultsFolder, "/all_threeprime.bed")

write.table( all.junctions, all.file, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t" )
write.table( all.threeprime, all.threeprime.file, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t" )
write.table( all.fiveprime, all.fiveprime.file, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t" )

print( "BedTools intersect junctions with list of known splice sites")
# first match junctions 
all.introns.cmd <- paste0("bedtools intersect -a ", all.file, " -b ", all_introns, " -wa -wb -loj -f 1" )
  
all.introns_intersect <- fread(all.introns.cmd)


# intersect with bedtools to find the annotations of each splice site
threeprime.cmd <- paste0( "bedtools intersect -a ", all.threeprime.file, " -b ",threeprime_file, " -wa -wb -loj -f 1" )

threeprime_intersect <- fread(threeprime.cmd)

fiveprime.cmd <- paste0( "bedtools intersect -a ", all.fiveprime.file, " -b ", fiveprime_file, " -wa -wb -loj -f 1" )

fiveprime_intersect <- fread(fiveprime.cmd)

# remove temporary files
rm.cmd <- paste("rm ", all.file, all.fiveprime.file, all.threeprime.file) 
system(rm.cmd)

# now I have two lists of splice site annotation
# for testing
#cluster <- all.introns[ all.introns$clusterID == "clu_4879" , ]

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
  cluster <- all.introns[ all.introns$clusterID == clu , ]
  
  # for each intron in the cluster:
  #   create vector of overlapping splice sites, indexed by the row of the intersect
  # five prime splice sites
  fprime <- apply( cluster, MAR = 1, FUN = function(x) {
    chr <- which( names(cluster) == "chr" )
    start <- which( names(cluster) == "start" )
    fiveprime_intersect[   
      fiveprime_intersect$V1 == x[chr] & 
      fiveprime_intersect$V2 == as.numeric( x[start] ),]
  } )
  # three prime splice sites
  tprime <- apply( cluster, MAR = 1, FUN = function(x) {
    chr <- which( names(cluster) == "chr" )
    end <- which( names(cluster) == "end" )
    threeprime_intersect[   
      threeprime_intersect$V1 == x[chr] & 
      threeprime_intersect$V2 == as.numeric( x[end] ),]
  } )

  # both splice sites
  bothSS <-  apply( cluster, MAR = 1, FUN = function(x) {
      chr <- which( names(cluster) == "chr" )
      start <- which(names(cluster) == "start")
      end <- which( names(cluster) == "end" )
      all.introns_intersect[   
        all.introns_intersect$V1 == x[chr] &
          all.introns_intersect$V2 == x[start] &
          all.introns_intersect$V3 == as.numeric( x[end] ) &
          all.introns_intersect$V4 == clu,]
    } )

  # find gene and ensemblID by the most represented gene among all the splice sites
  cluster_genes <- names(sort(table(do.call( what = rbind, tprime )$V8), decreasing = TRUE ))

  cluster_gene <- cluster_genes[ cluster_genes != "." ][1]
  # if no cluster gene found then leave as "."
  if( length(cluster_gene) == 0){
    cluster_gene == "."
  }
  # do the same for EnsemblID
  cluster_ensemblIDs <- names(sort(table(do.call( what = rbind, tprime )$V9), decreasing = TRUE ))
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
    coord[intron] <- paste(cluster[intron]$chr,cluster[intron]$start, cluster[intron]$end )

    gene[intron] <- cluster_gene
    ensemblID[intron] <- cluster_ensemblID

    # for each intron create vector of all transcripts that contain both splice sites
    transcripts[[intron]] <- unique( intersect( tprime[[intron]]$V10, fprime[[intron]]$V10 ) ) 

    verdict[intron] <- "error"
    if( # if neither are annotated
    all( tprime[[intron]]$V5 == ".") & all( fprime[[intron]]$V5 == "." )
    ){ verdict[intron] <- "cryptic_unanchored"
    }
    if( # if only one is annotated
    all( tprime[[intron]]$V5 == ".") & all( fprime[[intron]]$V5 != "." )
    ){ verdict[intron] <- "cryptic_threeprime"
    }
    if(
    all( tprime[[intron]]$V5 != ".") & all( fprime[[intron]]$V5 == "." )
    ){ verdict[intron] <- "cryptic_fiveprime"
    }
    if( # if both splice sites are annotated
      all( tprime[[intron]]$V5 != "." ) & all( fprime[[intron]]$V5 != "." )
    ){ 
      # test if the splice sites are paired in a known intron
      if( all(bothSS[[intron]]$V5 != ".") ){
        verdict[intron] <- "annotated"
      }else{ # both are annotated but never in the same junction
        verdict[intron] <- "skiptic"
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

    if( verdict.list[[clu]][2] == "skiptic" & verdict.list[[clu]][1] == "annotated" & verdict.list[[clu]][3] == "annotated"  ){
      classification.list[[clu]] <- paste0( classification.list[[clu]], " - skiptic")
    }

  }
  
}

print("Preparing results")

# match the lists together
all.introns$verdict <- unlist(verdict.list)[ match( paste( all.introns$chr, all.introns$start, all.introns$end ), unlist(coord.list)) ]

all.introns$gene <- unlist(gene.list)[ match( paste( all.introns$chr, all.introns$start, all.introns$end ), unlist(coord.list)) ]

all.introns$ensemblID <- unlist(ensemblID.list)[ match( paste( all.introns$chr, all.introns$start, all.introns$end ), unlist(coord.list)) ]

all.introns$transcripts <- unlist( transcripts.list )[ match( paste( all.introns$chr, all.introns$start, all.introns$end ), unlist(coord.list)) ]

all.introns$constitutive.score <-  unlist( constitutive.list )[ match( paste( all.introns$chr, all.introns$start, all.introns$end ), unlist(coord.list)) ]

#all.introns$prediction <-  unlist( classification.list )[ match( paste( all.introns$chr, all.introns$start, all.introns$end ), unlist(coord.list)) ]

# replace NA values with "."
all.introns$gene[ is.na( all.introns$gene) ] <- "."
all.introns$ensemblID[ is.na( all.introns$ensemblID) ] <- "."
# replace missing transcripts with "."
all.introns[ all.introns$transcripts == "", ]$transcripts <- "."
all.introns$constitutive.score <- signif(all.introns$constitutive.score, digits = 2)

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
  if( any(grepl( "cryptic", cluster$verdict)) | any( grepl("skiptic", cluster$verdict)) ){
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

# write out 
cluster_results <- paste0(resultsFolder,"/per_cluster_results.tab")
intron_results <- paste0(resultsFolder, "/per_intron_results.tab")

write.table( all.clusters, cluster_results, quote = FALSE, row.names = FALSE, sep = "\t" )
write.table( all.introns, intron_results, quote = FALSE, row.names = FALSE, sep = "\t" )

# prepare for PCA
counts <- counts[,meta$sample]
meta$group <- as.factor(meta$group)


make_pca <- function(counts,meta){
  dev <- apply( counts, MAR = 1, FUN = sd )
  # remove rows with 0 variance
  counts <- counts[ dev != 0, ]
  pca <- prcomp( t(counts), scale = TRUE )
  importance <- signif( summary(pca)$importance[2,], digits = 2) * 100
  pca <- as.data.frame(pca$x)
  #names(pca) <- paste0( names(pca), " (", signif( importance, digits =  2) * 100, "%)" )
  pca$groups <- meta$group[ match( rownames(pca), meta$sample )]
  return(list( pca, importance) )
}

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
                         nrow(all.introns[ all.introns$verdict == "skiptic",])  
                  )
                )
    return( summary )
}

# create all the objects for visualisation

clusters <- fix_clusters(all.clusters)
introns <- fix_introns(all.introns)
intron_summary <- intron_summary(all.introns)
cluster_summary <- cluster_summary(all.clusters) 
introns_to_plot <- leafcutter::get_intron_meta(rownames(counts))
cluster_ids <- introns_to_plot$clu 
pca <- make_pca(counts,meta)

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
      file = paste0( resultsFolder, "/",code,"_results.Rdata")
)


quit()


