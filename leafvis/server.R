library(shiny)
library(dplyr)
library(ggplot2)
library(DT)
library(leafcutter)
library(reshape2)
library(gridExtra)
library(intervals) # needed for pretty strand arrow placement
library(foreach)
library(shinycssloaders)
#data.table is required
# library(shinyjs)

filter_intron_table <- function(introns, clu, toSave=FALSE){
  d <- dplyr::filter(introns, clusterID == clu) %>% 
    dplyr::select( -clusterID, -gene, -ensemblID, -transcripts) %>% 
    arrange( desc(abs(deltapsi)))
  if( !toSave ){
    d <- rename(d, "ΔPSI" = deltapsi )
  }else{
    d <- rename(d, "dPSI" = deltapsi ) # fudge as grid arrange doesn't like greek letters
  }
  row.names(d) <- letters[1:nrow(d)] # letters is just a:z
  return(d)
}

# have a default option - this will be GTEX
# so the same server.R script can be given an Rdata object 
## by run_leafvis.R
# or can be run alone and it will load the default dataset.
# if( !exists("resultsData") ){
#   print("loading example dataset")
#   #defaultData <- "data/prudencio_fc_c9_results.Rdata"
#   defaultData <- "example/Brain_vs_Heart_results.Rdata"
#   if( file.exists(defaultData)){
#     load(defaultData)
#     # gencode_exons <- "data/gencode_hg38_all_exons.txt"
#     # exons_table <- as.data.frame(data.table::fread(gencode_exons))
#   }else{
#     stop("no dataset selected")
#   }
# }

#############
# SHINY APP
#############

server <- function(input, output) {
  output$logo <- renderImage({
    list( src = leafcutter_logo, alt = "", width = "15%", height = "15%" )
  }, deleteFile = FALSE)
  
  output$all_clusters <- DT::renderDataTable({
    datatable( clusters,
              escape = FALSE,
              rownames = FALSE,
              selection = 'single', 
              caption = "all significant clusters. N: number of introns within a cluster",
              fillContainer = FALSE,
              options = list(
                columnDefs = list(list(className = 'dt-center', targets = 0:5) )
              )
              ) 
  })
  output$sample_table <- DT::renderDataTable({
    datatable(sample_table,
              escape = FALSE,
              rownames = FALSE,
              fillContainer = FALSE,
              options <- list( searching = FALSE, paging = FALSE, info = FALSE )
              )
  })
  output$experimentCode <- renderText({
    paste("Experiment code:", code )
  })
  output$annotationCode <- renderText({
    paste("Annotation source:", basename(annotation_code) )
  })
  
  output$cluster_summary <- DT::renderDataTable({
    datatable(cluster_summary,
              escape = FALSE,
              rownames = FALSE,
              fillContainer = FALSE,
              options <- list( searching = FALSE, paging = FALSE, info = FALSE )
              )
  })
  
  output$intron_summary <- DT::renderDataTable({
    datatable(intron_summary,
              escape = FALSE,
              rownames = FALSE,
              fillContainer = FALSE,
              options <- list( searching = FALSE, paging = FALSE, info = FALSE )
              )
  })
  output$cluster_view = DT::renderDataTable({
    clu <- mycluster()
    if(!is.null(clu)){
      if(length(introns)){
        datatable( filter_intron_table(introns, clu, toSave=FALSE),
                 autoHideNavigation = TRUE, rownames = TRUE,
                 options <- list( searching = FALSE, paging = FALSE, info = FALSE)
          )
      }
    }
    
  })
  
  mygene <- eventReactive(input$all_clusters_rows_selected,{
    sel <- input$all_clusters_rows_selected
    if(is.null(sel)){ return(NULL)}
    gene  <- clusters[ sel, ]$gene
    gene <- gsub("<.*?>", "", gene) # strip out html italic tags
    return(gene)
  })
  
  mycluster <- eventReactive(input$all_clusters_rows_selected,{
    sel <- input$all_clusters_rows_selected
    if(is.null(sel)){return(NULL)}
    clusterID <- clusters[ sel, ]$clusterID
    return(clusterID)
  } )
   

  # NEW PLOTTING FUNCTIONS
  output$select_cluster_plot <- renderPlot({
    suppressWarnings( print(
      make_cluster_plot( mycluster(),
                       main_title = NA,
                       meta = meta,
                       cluster_ids = cluster_ids,
                       exons_table = exons_table,
                       counts = counts,
                       introns = introns)
    ))
  }, width = "auto", height = "auto",  res = 90
  )
  
  selectGenePlotInput <- function(all=FALSE){
    gene <- mygene()
    clu <- mycluster()
    if( ! is.null( gene) ){
      print( paste0( "gene is: ", gene))
      if(all != TRUE){
        clusterID <- clu
      }else{
        clusterID <- NULL
      }
      make_gene_plot(gene, counts = counts, introns = introns, exons_table = exons_table, cluster_list = clusters, clusterID = clusterID, introns_to_plot = introns_to_plot)
    }else{
      NULL
    }
  }
  
  selectClusterPlotInput <- function(title=NA){
    clu <- mycluster()
    if( !is.null(clu)) {
      make_cluster_plot( clu,
                       main_title = title,
                       meta = meta,
                       cluster_ids = cluster_ids,
                       exons_table = exons_table,
                       counts = counts,
                       introns = introns
                       )
    }else{
      NULL
    }
  }
  
  output$select_gene_plot <- renderPlot({
    suppressWarnings( print( selectGenePlotInput(all=FALSE) ) )
  }, width = "auto", height ="auto", res = 90
  )
  

  output$gene_title <- renderText({
    g <- mygene()
    if( is.null(mygene())){return("Cluster view\nSelect a cluster from the results table")}
    return( as.character( g )  ) 
  })

  output$cluster_title <- renderText({
    return( as.character(mycluster()))
  })
  
  # DOWNLOAD HANDLING
  output$downloadClusterPlot <- downloadHandler(
    filename = function() { paste0(mygene(),"_", mycluster(), '.pdf') },
    content = function(file) {
      plotTitle <- paste(mygene(), mycluster() )
      ggsave(file, plot = selectClusterPlotInput(title = plotTitle ), device = "pdf", width = 10, height = 5 )
    }
  )
  
  output$downloadClusterPlotWithTable <- downloadHandler(
    filename = function() { paste0(mygene(),"_", mycluster(), '_table.pdf') },
    content = function(file) {
      plotTitle <- paste(mygene(), mycluster() )
      clusterPlot <- selectClusterPlotInput(title=plotTitle)
      tablePlot <- tableGrob(filter_intron_table(introns, mycluster(), toSave=TRUE) )
      ggsave(file, plot = grid.arrange(clusterPlot, tablePlot, nrow =2),
           device = "pdf", width = 10, height = 8 )
    }
  )
  
  
  output$downloadGenePlot <- downloadHandler(
    filename = function() { paste( mygene(),"_","allClusters", '.pdf', sep='') },
    content = function(file) {
      ggsave(file, plot = selectGenePlotInput(all=TRUE), device = "pdf", width = 30, height = 5, limitsize = FALSE)
    }
  )
  # PCA
  
  output$pca_choices <- renderUI({
    choices <- names(pca[[1]])[ 1:(length(pca[[1]]) - 2) ]
    selectInput( inputId = "first_PC", label = "First principal component", choices = choices, selected = choices[1]  )
  })
  
  createPCAPlot <- function(){
    if( is.null(input$first_PC) ){
      return(NULL)
    }else{
    first_PC <- input$first_PC
    print(first_PC)
    second_PC <- names(pca[[1]])[ which( names(pca[[1]]) == first_PC) + 1 ]
    xlab <- paste0( first_PC, " (", pca[[2]][ which(names(pca[[1]]) == first_PC  ) ], "%)"  )
    ylab <- paste0( second_PC, " (", pca[[2]][ which(names(pca[[1]]) == first_PC ) + 1  ], "%)"  )
    
    pca_plot <- ggplot( pca[[1]], 
                        aes_string(y = first_PC,
                                  x = second_PC, 
                                  colour = "groups" ) ) + geom_point(size = 60 / nrow(pca[[1]])  ) +
      xlab( xlab ) +
      ylab( ylab )
    
    pca_plot
    }
  }
  
  output$pca_plot <- renderPlot({
    if( ! is.null( input$first_PC)){
      createPCAPlot()
    }else{
      NULL
    }
  }, width = "auto", height ="auto", res = 90)

  output$downloadPCAPlot <- downloadHandler(
    filename = function() { paste0('PCA.pdf') },
    content = function(file) {
    ggsave(file, plot = createPCAPlot(), device = "pdf", width = 7, height = 7 )
  })

}