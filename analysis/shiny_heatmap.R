library(shiny)
library(ComplexHeatmap)
library(circlize)

server <- function(input, output, session){
  
  # file upload
  filedata <- reactive({
    infile <- input$filename
    if (is.null(infile)){
      return(NULL)
    }
    read.csv(infile$datapath,sep = ",", header = T, row.names = 1)
  })
  
  # choose the max and min number in filedata
  observe({
    if (is.null(filedata())){
      return(NULL)
    }else{
      df <- filedata()
      min_value <- round(min(df))
      max_value <- round(max(df))
      updateNumericInput(session, "lowvalue", value = min_value)
      updateNumericInput(session, "highvalue", value = max_value)
    }
  })
  
  # choose the fit width and height for heatmap
  observe({
    if (input$downlist == "pdf"){
      updateNumericInput(session, "simpleheigh", value = 30)
      updateNumericInput(session, "simplewidth", value = 6)
    }else{
      updateNumericInput(session, "simpleheigh", value = 3000)
      updateNumericInput(session, "simplewidth", value = 4000)
    }
  })
  
  observe({
    if (input$complextype == "pdf"){
      updateNumericInput(session, "complexheigh", value = 30)
      updateNumericInput(session, "complexwidth", value = 6)
    }else{
      updateNumericInput(session, "complexheigh", value = 3000)
      updateNumericInput(session, "complexwidth", value = 4000)
    }
  })  
  
  ##########################################complex heatmap#######################################
  
  # prepare pre-data
  show_rowname <- reactive({
    rowname <- input$ori_rowname
    if (rowname == "none")
      return("FALSE")
    else
      return("TRUE")
  })
  
  show_columnname <- reactive({
    columnname <- input$ori_columnname
    if (columnname == "none")
      return("FALSE")
    else
      return("TRUE")
  })
  
  row_side <- reactive({
    #    rowside <- input$ori_rowname
    if (input$ori_rowname == "none")
      return(NULL)
    else
      input$ori_rowname
  })
  
  column_side <- reactive({
    columnside <- input$ori_columnname
    if (columnside == "none")
      return(NULL)
    else
      input$ori_columnname
  })
  
  # complex heatmap
  complexheatmap <- function(){
    if (!is.na(input$lowvalue) & !is.na(input$highvalue)){
      df <- filedata()
      df <- as.data.frame(df)
      Heatmap(df,
              name = input$nametxt,
              show_heatmap_legend = input$requirename,
              show_row_names = show_rowname(),
              show_column_names = show_columnname(),
              row_names_side = row_side(),
              column_names_side = column_side(),
              col = colorRamp2(c(input$lowvalue, 0, input$highvalue), c(input$lowcol, "white", input$highcol)),
              row_dend_width = unit(input$width_rowdend, "cm"),
              column_dend_height = unit(input$height_columndend, "cm"),
              clustering_distance_rows = input$distance,
              clustering_distance_columns = input$distance,
              clustering_method_rows = input$algorithm
      ) 
    }
  }
  
  output$complexmap <- renderPlot({
    if (!is.null(filedata()))
      complexheatmap()
  })
  
  output$complexgenelist <- renderTable({
    if (!is.null(filedata())){
      if (!is.na(input$lowvalue) & !is.na(input$highvalue)){
        px <- complexheatmap()
        dfx <- as.data.frame(filedata())
        listx <- row.names(dfx)[row_order(px)[[1]]]
        head(listx, 10)
      }
    }
  })
  
  # complex heatmap download
  output$downloadcomplex <- downloadHandler(
    filename <- function(){
      paste("cluster", input$complextype, sep = ".")
    },
    
    content <- function(file){
      if (input$complextype == "png")
        png(file, width = input$complexwidth, height = input$complexheigh, units ="px", res = 300)
      else
        pdf(file, width = input$complexwidth, height = input$complexheigh)
      print(complexheatmap())
      dev.off()
    }
    
  )  
  
  
  # complex gene list download
  output$downloadcomplexlist <- downloadHandler(
    filename <- "gene_list.csv",
    content <- function(file){
      if (!is.null(filedata())){
        px <- complexheatmap()
        dfx <- as.data.frame(filedata())
        listx <- row.names(dfx)[row_order(px)[[1]]]
        listx <- as.data.frame(listx)
        names(listx) <- "gene_list"
        write.csv(listx, file,row.names = FALSE,quote = FALSE)
      }
    }
  )  
  
}

ui <- dashboardPage(
  dashboardHeader(title = "ComplexHeatmap"),
  dashboardSidebar(
    fileInput("filename","Choose File to Upload:", accept = c(".csv")),
    sidebarMenu(
      menuItem("Quickly visualize", tabName = "quicklyplot"),
      menuItem("Complex visualize", tabName = "complexplot")
    )
  ),
  dashboardBody(
    tabItems(
      tabItem("quicklyplot",
              fluidRow(
                column(width = 9,
                       box(width = NULL,solidHeader = TRUE,
                           plotOutput("simpleplot", height = 500)
                       ),
                       fluidRow(
                         column(width = 6,
                                box(width = NULL, status = "info",
                                    title = "The gene list on the row",
                                    tableOutput("genelist") 
                                )
                         ),
                         column(width = 5, offset = 1,
                                box(width = NULL, status = "warning",
                                    downloadButton("downloadlist",label = "Download Gene List")    
                                )
                         )
                       )
                ),
                column(width = 3,
                       box(width = NULL,status = "warning",
                           numericInput("simpleheigh",label = "Graph heigh value",value = 3000),
                           numericInput("simplewidth",label = "Graph width value",value = 4000),
                           radioButtons("downlist",label = "Select the graph type",,choices = c("png","pdf")),
                           downloadButton("downloadsimple",label = "Download Heatmap")
                       )
                )
              )
      ),
      tabItem("complexplot",
              fluidRow(
                column(width = 8,
                       box(width = NULL,solidHeader = TRUE,
                           plotOutput("complexmap",height = 500)
                       ),
                       fluidRow(
                         column(width = 6,
                                box(width = NULL, status = "info",
                                    title = "The gene list on the row",
                                    tableOutput("complexgenelist")
                                )
                         ),
                         column(width = 5, offset = 1,
                                box(title = "Download",
                                    solidHeader = T, status = "info",
                                    width = NULL,
                                    h5("Choose height and width for heatmap"),
                                    fluidRow(
                                      column(width = 6,
                                             numericInput("complexheigh",label = "Graph heigh value",value = 3000)
                                      ),
                                      column(width = 6,
                                             numericInput("complexwidth",label = "Graph width value",value = 4000)
                                      ),
                                      column(width = 12,
                                             radioButtons("complextype",label = "Select the graph type",choices = c("png","pdf"))
                                      )
                                    )
                                ),
                                
                                box(width = NULL,status = "warning",
                                    downloadButton("downloadcomplexlist",label = "Download Gene List"),
                                    downloadButton("downloadcomplex",label = "Download Heatmap")
                                )
                         )
                       )
                ),
                
                
                box(title = "The parameter for complex heatmap",
                    solidHeader = T, status = "info",
                    collapsible = T, collapsed = F,
                    width = 4,
                    
                    fluidRow(
                      box(title = "Expression To Color Option",
                          solidHeader = T, status = "info",
                          width = 12,
                          h5("The left is expression value, the right is color"),
                          fluidRow(
                            column(width = 6,
                                   numericInput("lowvalue", label = "Low Value", value = NULL,-5,0)
                            ),
                            column(width = 6,
                                   selectInput("lowcol", label = "Low color", choices = c("green","blue", "purple", "red", 
                                                                                          "orange", "yellow", "white"),
                                               selected = "green")
                            ),
                            column(width = 6,
                                   numericInput("highvalue", label = "High Value", value = NULL,0,5)
                            ),
                            column(width = 6,
                                   selectInput("highcol", label = "High color", choices = c("red", "orange", "yellow", 
                                                                                            "green", "blue", "purple", "white"),
                                               selected = "red")
                            )
                          )
                      ),
                      
                      box(title = "Heatmap Name",
                          solidHeader = T, status = "info",
                          width = 12,
                          h5("Please choose the name option for heatmap"),
                          fluidRow(
                            column(width = 6,
                                   selectInput("requirename", label = "Legend or not", c("TRUE", "FALSE"))
                            ),
                            column(width = 6,
                                   textInput("nametxt", label = "Legend Name", value = "expression")
                            )
                          )
                      ),
                      
                      box(title = "Row Relative Option",
                          solidHeader = T, status = "info",
                          width = 12,
                          h5("Choose orientation for row name and the height for dend"),
                          fluidRow(
                            column(width = 6,
                                   selectInput("ori_rowname", label = "Select orientation for row name",
                                               c("right","left","none"))
                            ),
                            column(width = 6,
                                   numericInput("width_rowdend",label = "Select height for row dend",
                                                value = 4,1,10)
                            )
                          )
                      ),
                      
                      box(title = "column Relative Option",
                          solidHeader = T, status = "info",
                          width = 12,
                          h5("Choose orientation for column name and the height for dend"),
                          fluidRow(
                            column(width = 6,
                                   selectInput("ori_columnname", label = "Select orientation for column name",
                                               c("bottom","top","none"))
                            ),
                            column(width = 6,
                                   numericInput("height_columndend",label = "Select height for column dend",
                                                value = 2,1,10)
                            )
                          )                          
                      ),
                      
                      box(title = "Cluster Relative Option",
                          solidHeader = T, status = "info",
                          width = 12,
                          h5("Choose distance and algorithm to cluster"),
                          fluidRow(
                            column(width = 6,
                                   selectInput("distance", label = "Select distance to cluster",
                                               c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski"))
                            ),
                            column(width = 6,
                                   selectInput("algorithm", label = "Select algorithm to cluster",
                                               c("complete", "single", "average", "centroid", "median"))
                            )
                          )
                      )
                    )
                )
              )
      )
    )
  )
)

shinyApp(ui, server)

  