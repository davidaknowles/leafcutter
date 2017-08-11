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

if (!exists("introns")) load("example/Brain_vs_Heart_results.Rdata")

#############
# SHINY APP
#############

server <- function(input, output) {
  output$logo <- renderImage({
    list( src = leafcutter_logo, alt = "", width = "15%", height = "15%" )
  }, deleteFile = FALSE)
  
  output$all_clusters <- DT::renderDataTable({
    datatable( clusters[,c("gene","coord","N","FDR","annotation")],
              escape = FALSE,
              rownames = FALSE,
              colnames = c('Genomic location'='coord','Gene'='gene','N'='N','Annotation'='annotation','q'='FDR'),
              selection = 'single', 
              caption = "Click on a row to plot the corresponding visualization. N: number of introns within a cluster. q: Benjamini–Hochberg q-value.",
              fillContainer = FALSE,
              options = list(
                columnDefs = list(list(className = 'dt-center', targets = 2:4) )
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
  
  mycoord <- eventReactive(input$all_clusters_rows_selected,{
    sel <- input$all_clusters_rows_selected
    if(is.null(sel)){return("")}
    coord <- clusters[ sel, ]$coord
    return(coord)
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
  
  # UCSC LINKS
  output$viewClusterUCSC <- renderUI({
    db <- strsplit(basename(annotation_code), split = "_")[[1]][2]
    coord <- mycoord()
    # zoom out the position a bit
    chr <- strsplit(coord, ":")[[1]][1]
    start <- as.numeric(strsplit( strsplit(coord, ":")[[1]][2], "-" )[[1]][1])
    end <- as.numeric(strsplit( strsplit(coord, ":")[[1]][2], "-" )[[1]][2])
    start <- start - 100
    end <- end + 100
    coord <- paste0(chr, ":", as.character(start), "-", as.character(end))
    url <- paste0( "http://genome.ucsc.edu/cgi-bin/hgTracks?&db=",db,"&position=", coord )
    return(tags$a(href = url, "view on UCSC", target = "_blank", class = "btn btn_default shiny-download-link  shiny-bound-output", id = "UCSC" ) )
    })
  
  output$viewGeneUCSC <- renderUI({
    db <- strsplit(basename(annotation_code), split = "_")[[1]][2]
    gene <- mygene()
    url <- paste0( "http://genome.ucsc.edu/cgi-bin/hgTracks?&db=",db,"&singleSearch=knownCanonical&position=", gene)
    return(tags$a(href = url, "view on UCSC", target = "_blank", class = "btn btn_default", id = "UCSC" ) )
  })
  
  # PCA
  
  output$pca_choices <- renderUI({
    choices <- names(pca[[1]])[ 1:(length(pca[[1]]) - 2) ]
    selectInput( inputId = "first_PC", label = "First principal component", choices = choices, selected = choices[1]  )
  })
  
  createPCAPlot <- function(){
    if( is.null(input$first_PC) ){
      return(NULL)
    }else{
    first_PC <- input$first_P
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