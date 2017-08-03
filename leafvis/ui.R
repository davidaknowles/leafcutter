library(shiny)
library(DT)
#library(shinycssloaders)

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
              h4(id = "title","cluster results"),
              hr(),
              div(
                # withSpinner(DT::dataTableOutput("all_clusters"))
                DT::dataTableOutput("all_clusters")
              )
            )
          ),
        column(6,
          
          ### CLUSTER VIEW
          
          div(id="clusterView",
            h4(id="title","cluster view"),
            hr(),
            h4(id="title", strong(  em( textOutput("gene_title") ) ), textOutput("cluster_title"), align = "left"),
            div(
              div(id = "welcome", 
                    h3("Welcome to the Leafcutter visualisation app!"),
                    p(HTML(paste0("To visualise a cluster, click a row in ",strong("cluster results.") ) ) ),
                    p(HTML(paste0("All clusters found within a gene are visualised in the ", strong("gene view"), " below.") ) )
                  ),
              #withSpinner(plotOutput("select_cluster_plot", width = "100%") )
              plotOutput("select_cluster_plot", width = "100%")
            ),
            DT::dataTableOutput("cluster_view"),
              hr(),
              div(class = "download_btn",
                downloadButton("downloadClusterPlot", label = "save plot", class = NULL),
                downloadButton("downloadClusterPlotWithTable", label = "save plot + table", class = NULL)
              )
          )
        )
      ),
      
      ### GENE VIEW
      
      hr(),
      verticalLayout(
        div(id="geneView",
          h4(id="title","gene view"),
          hr(),
          # div(id = "welcome", 
          #   h3("No gene selected")
          #   ),
          #withSpinner(plotOutput("select_gene_plot", width="100%")),
          plotOutput("select_gene_plot", width="100%"),
          hr(),
          div(class = "download_btn",
            downloadButton("downloadGenePlot", label = "Save plot", class = NULL)
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
    # tabPanel("Settings",
    #   fluidRow()
    # ),
    collapsible = TRUE,
    inverse = FALSE,
    position = "fixed-top",
    footer = div(id = "footer",
      p("created by Jack Humphrey", "(https://github.com/jackhump)"),
      tags$head(tags$link(rel="shortcut icon", href="favicon.ico"))      
    )

  )
)



