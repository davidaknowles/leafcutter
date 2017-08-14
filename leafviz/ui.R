library(shiny)
library(DT)
library(shinycssloaders)

ui <- shinyUI(
  navbarPage(title = "Leafcutter",
    tabPanel("All clusters",
      tags$style(type="text/css", "
body {
    padding-top: 70px;
}
#title {
  color: black;
  margin-top:5px;
}
hr {
    margin-top: 10px;
    margin-bottom: 10px;
}
#geneView {
  border-radius: 25px;
  border-color: whitesmoke;
  border-style: solid;
  border-width:5px;
  padding: 10px;
}
.navbar {
background-size: 50px 80%;
background-image: url(leafcutter_app_logo.png);
background-repeat: no-repeat;
background-position: 2%;
padding-left: 2%;
}

.download_btn{
  margin-top: 10px;
  margin-bottom: 10px;
  display: block;
  text-align: center;
}

#clusterView {
  padding-top: 10px;
  padding-left: 5px;
  padding-right: 5px;
  border-radius: 25px;
  border-color: whitesmoke;
  border-style: solid;
  border-width:5px;
}

#title {
  text-align: center;
}

#welcome {
  text-align: center;
  font-size: large;
  position: absolute;
  padding: 15%;
}

#summary {
  text-align: center;
}

#footer {
  margin: auto;
  padding: 10px;
  text-align: center;
}

#UCSC {
  border-color: #ccc;
  color: #333;
}

#UCSC:hover {
    background-color: #E6E6E6;
}

"),
        fluidRow(
          div(
            column(6, 
                  style="background-color: WhiteSmoke;
                         border-radius: 25px;
                         padding: 10px;
                      
                         margin-top: 0px;
                        ",
              h4(id = "title","Differential splicing events (clusters)"),
              hr(),
              div(
                withSpinner(DT::dataTableOutput("all_clusters"))
              )
            )
          ),
        column(6, 
          
          ### CLUSTER VIEW
          
          div(id="clusterView",
            h4(id="title","Splicing event visualization"),
            hr(),
            h4(id="title", strong(  em( textOutput("gene_title") ) ), textOutput("cluster_title"), align = "left"),
            div(
              div(id = "welcome", 
                    h3("Welcome to the Leafcutter visualisation app!"),
                    p(HTML(paste0("To visualise a cluster, click a row in ",strong("cluster results.") ) ) ),
                    p(HTML(paste0("All clusters found within a gene are visualised in the ", strong("gene view"), " below.") ) )
                  ),
              withSpinner(plotOutput("select_cluster_plot", width = "100%") )
            ),
            DT::dataTableOutput("cluster_view"),
              hr(),
              div(class = "download_btn",
                downloadButton("downloadClusterPlot", label = "save plot", class = NULL),
                downloadButton("downloadClusterPlotWithTable", label = "save plot + table", class = NULL),
                htmlOutput("viewClusterUCSC", inline = TRUE)
              )
          )
        )
      ),
      
      ### GENE VIEW
      
      hr(),
      verticalLayout(
        div(id="geneView",
          h4(id="title","Gene-level visualization"),
          hr(),
          # div(id = "welcome", 
          #   h3("No gene selected")
          #   ),
          withSpinner(plotOutput("select_gene_plot", width="100%")),
          hr(),
          div(class = "download_btn",
            downloadButton("downloadGenePlot", label = "Save plot", class = NULL),
            htmlOutput("viewGeneUCSC", inline = TRUE)
          )
        )
      )
    ),
      
    tabPanel("Summary",
      fluidRow(
        column(6,offset=3,
          # this is where a summary table goes counting all the significant clusters and introns
          h3("Experiment",id="summary"),
          hr(),
          h5(textOutput("experimentCode"),id="summary"),
          h5(textOutput("annotationCode"),id="summary"),
          br(),
          h4("Samples",id="summary"),
          DT::dataTableOutput("sample_table"),
          br(),
          h4("Clusters",id="summary"),
          hr(),
          DT::dataTableOutput("cluster_summary"),
          br(),
          h4("Junctions",id="summary"),
          hr(),
          DT::dataTableOutput("intron_summary"),
          br(),
          br()
        )
      )
    ),
    tabPanel("PCA",
      fluidRow(
        column(6, offset =3,
          # plot different principal components of the splice junction counts
          br(),
          div(id = "download_btn", 
            uiOutput("pca_choices") 
            ),
          plotOutput("pca_plot")
        )
      ),
      div(class = "download_btn",
          downloadButton("downloadPCAPlot", label = "save plot", class = NULL)
      )
    ),
    # SETTINGS - any needed?
    tabPanel("About",
       fluidRow(
         column(
           6,
           offset=3,
           div(
            h2("What is this?"),
            p( "This R", tags$a(href="https://shiny.rstudio.com","Shiny", target = "_blank"), 
               "app presents and visualises the results of running",
               a(href="https://github.com/davidaknowles/leafcutter", "Leafcutter,", target = "_blank"),
               "a software tool that quantifies RNA-seq splicing in an annotation-free way.", 
               "Full documentation of the package is available", a(href="http://davidaknowles.github.io/leafcutter/", "here.", target = "_blank") 
              ),
            h2( "Differential splicing events"),
            p( "A cluster is defined as set of overlapping spliced junctions or introns.", 
               "Clusters are initially ranked in the cluster results table by adjusted P value."
               ),
            tags$ul(
              tags$li( strong("Gene -"), "the HUGO gene sympbol for that gene."),
              tags$li( strong("Genomic location - "), "the coordinate span of the largest intron in the cluster."  ),
               # tags$li( strong("clusterID - "), "the unique id assigned to the cluster by leafcutter." ), 
                tags$li( strong("N -") , "the number of introns in the cluster." ), 
                tags$li( strong("q -"), "the Benjamini-Hochberg adjusted P value of the multinomial test of intron counts between conditions."),
              tags$li( strong("Annotation -"), "whether every intron in the cluster is supported by annotation (annotated) or contains at least one unannotated junction (cryptic).")
                ),
            h2("Splicing event visualization"),
            p("To view a cluster, click on a row in the cluster results table. This will start the plotting function.", 
              "A cluster plot and table are generated each time a row in the cluster results is clicked."),
            p("For a chosen cluster, the mean number of splice junctions supporting each intron is calculated for both conditions and then normalised as a fraction of the total counts.",
              "Therefore for each condition the normalised counts will add up to 1.",
              "Each intron is plotted as a line connecting its start and end coordinates with a thickness proportional to the displayed normalised count value.",
              "The colour of the intron line indicates whether it is present in the annotation (red) or not (pink).",
              "Any exons that are annotated as flanking or being contained within the cluster are added as rectangles to the plot.",
              "If exons from multiple genes are connected by introns then their exons will be coloured according to their gene of origin."
              ),
            p("Each intron is presented as a row in the cluster view table."
              ),
           tags$ul(
             tags$li( strong("chr, start, end"), "the genomic coordinates of the intron." ), 
             tags$li( strong("verdict") , "the support given to two splice sites of the intron (start and end) by annotation." ),
             tags$ul(
                tags$li( em("annotated -"), "both splice sites are present in an annotated junction"),
                tags$li( em("novel annotated pair -"), "both splice sites are annotated but are not annotated as being paired in a junction"),
                tags$li( em("cryptic_fiveprime -"), "the 3\' splice site is annotated but the 5\' is not."),
                tags$li( em("cryptic_threeprime -"), "the 5\' splice site is annotated but the 3\' is not ")
             ),
             tags$li( strong("dPSI"), "The Leafcutter algorithm estimates a dPSI (delta percent spliced in) value for each intron and this is displayed in the cluster view table."  )
           ),
           p("Each intron in the cluster is ranked by the absolute dPSI value. The ranked introns are then assigned a letter value for labelling."),
           
           h2("Gene-level visualization"),
           p("This visualises all clusters discovered by Leafcutter that can be assigned to a particular gene.",
             "Exons are taken from the provided annotation and plotted as black rectangles.",
             "Each junction in each cluster is plotted as curved line with uniform thickness",
             "Junctions from significant clusters are coloured according to the estimated dPSI (see above), whereas junctions from clusters that are not significant are coloured grey.",
             "Note that the genomic coordinates are deliberately warped to give more space to the clusters."
             ),
            
            
              h3(tags$a(href="http://davidaknowles.github.io/leafcutter/articles/Visualization.html",
                     "How to visualise your own Leafcutter results", target = "_blank"))
            )
         )
        )
     ),
    collapsible = TRUE,
    inverse = FALSE,
    position = "fixed-top",
    footer = div(id = "footer",
      p("written by Jack Humphrey, David Knowles & Yang Li.",
        tags$a(href="https://github.com/davidaknowles/leafcutter/tree/master/leafvis", "Fork on GitHub.", target = "_blank")
    ),
      tags$head(tags$link(rel="shortcut icon", href="favicon.ico"))      
    )

  )
)



