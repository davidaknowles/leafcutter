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
library(grid)
library(gtable)
library(ggrepel)
#data.table is required
# library(shinyjs)
#### options for debugging
#options(shiny.trace=TRUE)
# options(shiny.reactlog=TRUE)
#load("example/Brain_vs_Heart_results.Rdata")
#source("../leafcutter/R/make_gene_plot.R")
#make_gene_plot("MICAL3", counts = counts, introns = introns, exons_table = exons_table, cluster_list = clusters, clusterID = "clu_36585", introns_to_plot = introns_to_plot)
#make_gene_plot("MICAL3",counts = counts, introns = introns, exons_table = exons_table, cluster_list = clusters, clusterID = NULL, introns_to_plot = introns_to_plot)
#source("../leafcutter/R/make_cluster_plot.R")
#make_cluster_plot( "clu_8845",main_title = c("RILP1", "clu_8845"), meta = meta, cluster_ids = cluster_ids, exons_table = exons_table, counts = counts,introns = introns)


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

getGeneLength <- function(gene_name){
  # gets length of gene in nucleotides and decides on a pixel length for the gene plot
  # RBFOX1 is 1.7Mbp - scale to 5000px
  # most genes are < 100kb
  exons <- exons_table[ exons_table$gene_name == gene_name, ]
  geneStart <- min(exons$start)
  geneEnd <- max(exons$end)
  geneLength <- geneEnd - geneStart
  #print(geneLength)
  if( geneLength >1E6){
    pixels <- 5000 # scales RBFOX1 to 5000px
  }
  if( geneLength > 5e5 & geneLength < 1e6){
    pixels <- 3000
  }
  if( geneLength > 1.5e5 & geneLength <= 5e5){
    pixels <- 2000
  }
  if( geneLength <= 1.5e5){
    pixels <- "auto"
  }
  #print(pixels)
  return(pixels)
}
# test
#getGeneLength("RBFOX1", 530)


if (!exists("introns")){
  load("example/Brain_vs_Heart_results.Rdata")
  defaultValue <- 12 #RBFOX1
  # for testing - simulate data aligned to genome without "chr" in chr name
  #introns_to_plot$chr <- gsub("chr","", introns_to_plot$chr)
  #row.names(counts) <- gsub("chr", "", row.names(counts))

}else{
  defaultValue <- NULL
}

#############
# SHINY APP
#############

server <- function(input, output, session) {
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
                pageLength = 15,
                columnDefs = list(list(className = 'dt-center', targets = 0:4) )
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

  # observeEvent(input$toggleInstruct, {
  #   renderUI()
  # })

  onclick("welcome", toggle(id = "popupInstruct", anim = TRUE) )

  observeEvent( input$aboutLink, {
    updateTabsetPanel(session, "navBarPage", selected = "About")
  })

  observeEvent( input$aboutLink2, {
    updateTabsetPanel(session, "navBarPage", selected = "About")
  })

  # SUMMARY - does this need to be in server and not just UI?

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

  # TABLES

  output$cluster_view = DT::renderDataTable({
    clu <- mydata()$cluster
    if(!is.null(clu)){
      if(length(introns)){
        datatable( filter_intron_table(introns, clu, toSave=FALSE),
                 autoHideNavigation = TRUE, rownames = TRUE,
                 options <- list( searching = FALSE, paging = FALSE, info = FALSE)
          )
      }
    }else{
      print("no cluster selected!")
    }

  })

  # SET REACTIVE VALUE WITH A DEFAULT

  values <- reactiveValues(default = defaultValue) # RBFOX1 in the Brain vs Heart dataset
  # REACTIVE VALUE IS UPDATED BY INPUT
  observeEvent(input$all_clusters_rows_selected,{
    #print("new row selected!")
    values$default <- input$all_clusters_rows_selected # if all_clusters_rows_selected changes then update value - this sets everything!
    #print(paste0("VALUE: ", values$default ))
  })

  # USE REACTIVE VALUE TO GENERATE ALL VARIABLES NEEDED

  mydata <- eventReactive(values$default,{
    sel <- values$default
    gene  <- clusters[ sel, ]$gene
    gene <- gsub("<.*?>", "", gene) # strip out html italic tags
    width <- getGeneLength(gene)
    clusterID <- clusters[ sel, ]$clusterID
    coord <- clusters[ sel, ]$coord
    return(list(gene = gene, width = width, cluster = clusterID, coord = coord) )
  })


  # PLOTTING

  output$select_cluster_plot <- renderPlot({
    plotTitle <- c(mydata()$gene, as.character(mydata()$cluster) )
    suppressWarnings( print(
      make_cluster_plot( mydata()$cluster,
                       main_title = plotTitle,
                       meta = meta,
                       cluster_ids = cluster_ids,
                       exons_table = exons_table,
                       counts = counts,
                       introns = introns)
    ))
  }, width = "auto", height = "auto",  res = 90
  )

  observeEvent(values$default,{
    output$select_gene_plot <- renderPlot({
    suppressWarnings( print(
      make_gene_plot(mydata()$gene, counts = counts, introns = introns, exons_table = exons_table, cluster_list = clusters, clusterID = mydata()$clusterID, introns_to_plot = introns_to_plot, debug=F)
      )
    )
    }, width = mydata()$width, height = "auto", res = 90 # try changing height param
   )
  })

  # TITLES

  output$gene_title <- renderText({
    return( as.character( mydata()$gene )  )
  })

  output$cluster_title <- renderText({
    return( as.character(mydata()$cluster))
  })

  # DOWNLOAD HANDLING

  output$downloadClusterPlot <- downloadHandler(
    filename = function() { paste0(mydata()$gene,"_", mydata()$cluster, '.pdf') },
    content = function(file) {
      plotTitle <- c(mydata()$gene, as.character(mydata()$cluster) )
      ggsave(file,
             plot = make_cluster_plot( mydata()$cluster,
                                  main_title = plotTitle,
                                  meta = meta,
                                  cluster_ids = cluster_ids,
                                  exons_table = exons_table,
                                  counts = counts,
                                  introns = introns),
             device = "pdf", width = 10, height = 5 )
    }
  )

  output$downloadClusterPlotWithTable <- downloadHandler(
    filename = function() { paste0(mydata()$gene,"_", mydata()$cluster, '_table.pdf') },
    content = function(file) {
      plotTitle <- c(mydata()$gene, as.character(mydata()$cluster ) )
      clusterPlot <- make_cluster_plot( mydata()$cluster,
                                        main_title = plotTitle,
                                        meta = meta,
                                        cluster_ids = cluster_ids,
                                        exons_table = exons_table,
                                        counts = counts,
                                        introns = introns)
      # make table theme
      tableTheme <- ttheme_minimal(
        core=list(bg_params = list(fill = c("whitesmoke","white"), col=NA)
        ),
        colhead=list(fg_params=list(col="black", fontface="bold"),
                     bg_params = list(fill="white")),
        rowhead=list(fg_params=list(col="black"),
                     bg_params = list(fill=c("white", "whitesmoke"))))

      mytable <- tableGrob(filter_intron_table(introns, mydata()$cluster, toSave=TRUE), theme = tableTheme )
      mycols <- ncol(mytable)
      mytable$widths <- unit( c( 1/(3*mycols), rep(1/mycols, mycols-1) ), "npc")

      mytable <- gtable_add_grob(mytable,
                           grobs = segmentsGrob( # line across the bottom
                             x0 = unit(0,"npc"),
                             y0 = unit(0,"npc"),
                             x1 = unit(1,"npc"),
                             y1 = unit(0,"npc"),
                             gp = gpar(lwd = 2.0)),
                           t = 2, b = nrow(mytable), l = 1, r = mycols)

      mytable <- gtable_add_grob(mytable,
                                 grobs = segmentsGrob( # line across the bottom
                                   x0 = unit(0,"npc"),
                                   y0 = unit(0,"npc"),
                                   x1 = unit(1,"npc"),
                                   y1 = unit(0,"npc"),
                                   gp = gpar(lwd = 2.0)),
                                 t = 1, b = 1, l = 1, r = mycols)

      ggsave(file, plot = grid.arrange(clusterPlot, mytable, nrow =2),
           device = "pdf", width = 10, height = 8 )
    }
  )

  output$downloadGenePlot <- downloadHandler(
    filename = function() { paste( mydata()$gene,"_","allClusters", '.pdf', sep='') },
    content = function(file) {
      ggsave(file,
             plot = make_gene_plot(mydata()$gene, counts = counts, introns = introns, exons_table = exons_table, cluster_list = clusters, clusterID = NULL, introns_to_plot = introns_to_plot),
             device = "pdf", width = ifelse( mydata()$width == "auto", yes = 10, no = mydata()$width / 100 ), height = 6, units = "in", limitsize = FALSE)
    }
  )

  # UCSC LINKS
  output$viewClusterUCSC <- renderUI({
    coord <- mydata()$coord
    db <- strsplit(basename(annotation_code), split = "_")[[1]][2]
    # guess species from genome build - if not possible then leave blank.
    org <- NULL
    if( grepl("hg", db)){ org <- "human"}
    if( grepl("mm", db)){ org <- "mouse"}
    if( is.null(org)){
      orgChoice <- ""
    }else{
      orgChoice <- paste0("&org=",org)
    }

    # zoom out the position a bit
    chr <- strsplit(coord, ":")[[1]][1]
    start <- as.numeric(strsplit( strsplit(coord, ":")[[1]][2], "-" )[[1]][1])
    end <- as.numeric(strsplit( strsplit(coord, ":")[[1]][2], "-" )[[1]][2])
    start <- start - 100
    end <- end + 100
    coord <- paste0(chr, ":", as.character(start), "-", as.character(end))
    url <- paste0( "http://genome.ucsc.edu/cgi-bin/hgTracks?", orgChoice, "&db=",db,"&position=", coord )
    return(tags$a(href = url, "view on UCSC", target = "_blank", class = "btn btn_default", id = "UCSC" ) )
    })

  output$viewGeneUCSC <- renderUI({
    db <- strsplit(basename(annotation_code), split = "_")[[1]][2]
    org <- NULL
    if( grepl("hg", db)){ org <- "human"}
    if( grepl("mm", db)){ org <- "mouse"}
    if( is.null(org)){
      orgChoice <- ""
    }else{
      orgChoice <- paste0("&org=",org)
    }
    gene <- mydata()$gene
    url <- paste0( "http://genome.ucsc.edu/cgi-bin/hgTracks?",orgChoice,"&db=",db,"&singleSearch=knownCanonical&position=", gene)
    return(tags$a(href = url, "view on UCSC", target = "_blank", class = "btn btn_default", id = "UCSC" ) )
  })

  #### PCA
  #names(pca[[1]]) <- gsub(" ", "_", names(pca[[1]]))
  output$pca_choices <- renderUI({
    choices <- names(pca[[1]])[ grepl("^PC[0-9]", names(pca[[1]])) ] # find all PCs
    selectInput( inputId = "first_PC", label = "First principal component", choices = choices, selected = choices[1]  )
  })
  output$pca_colour_choices <- renderUI({
    choices <- names(pca[[1]])[ !grepl("^PC[0-9]", names(pca[[1]])) ]
    selectInput( inputId = "colour_choice", label = "Colour points by", choices = choices, selected = choices[1]  )
  })
  output$pca_shape_choices <- renderUI({
    choices <- names(pca[[1]])[ !grepl("^PC[0-9]", names(pca[[1]])) ]
    selectInput( inputId = "shape_choice", label = "Shape points by", choices = choices, selected = choices[1]  )
  })


  createPCAPlot <- function(){
    if( is.null(input$first_PC) ){
      return(NULL)
    }else{
    first_PC <- input$first_PC
    colour_choice <- input$colour_choice
    shape_choice <- input$shape_choice
    #print(first_PC)
    second_PC <- names(pca[[1]])[ which( names(pca[[1]]) == first_PC) + 1 ]
    xlab <- paste0( first_PC, " (", pca[[2]][ which(names(pca[[1]]) == first_PC  ) ], "%)"  )
    ylab <- paste0( second_PC, " (", pca[[2]][ which(names(pca[[1]]) == first_PC ) + 1  ], "%)"  )

    pca_plot <- ggplot( pca[[1]],
                        aes_string(x = first_PC,
                                  y = second_PC,
                                  colour = colour_choice,
                                  shape = shape_choice) ) + geom_point(size = 60 / nrow(pca[[1]])  ) +
      xlab( xlab ) +
      ylab( ylab ) +
      theme_classic()

    pca_plot
    }
  }

  output$pca_plot <- renderPlot({
    if( ! is.null( input$first_PC)){
      createPCAPlot()
    }
  }, width = "auto", height ="auto", res = 90)

  output$downloadPCAPlot <- downloadHandler(
    filename = function() { paste0('PCA.pdf') },
    content = function(file) {
    ggsave(file, plot = createPCAPlot(), device = "pdf", width = 7, height = 7 )
  })

}
