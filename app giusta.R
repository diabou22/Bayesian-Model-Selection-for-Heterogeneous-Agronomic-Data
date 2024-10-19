library(shiny)
library(maps)
library(mapproj)
library(ggplot2)
library(rstan)
library(brms)
library(tidyverse)
library(bayesplot)
library(tidybayes)
library(dplyr)
library(readr)
library(DT)
library(shinythemes)
library(htmltools)
library(caret)

# Load data
load("file_dati_tesi.RData")

# Data for the table
trial_data <- data.frame(
  Site = c("Debrecen", "Debrecen", "Debrecen", "Gaillac", "Gaillac", "Gaillac", "Karlsruhe", "Karlsruhe", "Karlsruhe"),
  Variety_ID = c("B73", "F1808", "F7028", "B73", "F1808", "F7028", "B73", "F1808", "F7028"),
  Treatment_Rainfed = c(15, 9, 9, 12, 6, 6, 16, 8, 8),
  Treatment_Watered = c(6, 4, 4, 8, 4, 4, 12, 6, 6)
)

# Load models into the global environment
load("allmodels.RData")

ui <- fluidPage(
  theme = shinytheme("sandstone"),
  titlePanel("Bayesian Model Selection for Heterogeneous Agronomic Data"),
  navbarPage("by Diabou Mbaye",
             id = "navBar",
             
             # Data panel
             tabPanel("Data summary", value = "df_complete",
                      strong("Statistical summary of filtered dataset is shown below:"),
                      verbatimTextOutput("dataSummary")
             ),
             
             # Download panel
             tabPanel("Download", value = "df_tesi",
                      radioButtons(
                        inputId = "filetype",
                        label = "Select filetype:",
                        choices = c("csv", "tsv","txt", "docx"),
                        selected = "csv"
                      ),
                      checkboxGroupInput(
                        inputId = "selected_var",
                        label = "Select variables:",
                        choices = names(df_tesi),
                        selected = c("PL", "TH", "ear.height", "anthesis", "silking", 
                                     "anthesis.silking.interval", "PL", 
                                     "grain.number", "grain.yield", "TH",
                                     "grain.weight", "Site", "year", "Experiment", "plotID", "treatment", "Replicate",
                                     "block", "Row", "Column", "Accession", "Code_ID", "Variety_ID")
                      ),
                      dataTableOutput("filetable"),
                      downloadButton("download_data", "Download data")
             ),
             
             # Location panel
             tabPanel("Geolocation Map", value = "map",
                      sidebarLayout(
                        sidebarPanel(
                          h4("Filtered Dataset Table"),
                          tableOutput("trialTable")
                        ),
                        mainPanel(
                          tags$iframe(
                            src = "https://urgi.versailles.inra.fr/ephesis/ephesis/viewer.do#dataResults/trialSetIds=42",
                            width = "100%",
                            height = "600px"
                          )
                        ))
             ),
             
             # Distribution panel
             tabPanel("Distribution", value = "graph",
                      radioButtons("dist", "Distribution type:",
                                   c("Normal" = "norm",
                                     "Uniform" = "unif",
                                     "Log-normal" = "lnorm",
                                     "Exponential" = "exp")),
                      selectInput(inputId = "var", 
                                  label = "Choose a variable to display",
                                  choices = c("PL", "TH", "ear.height", "anthesis", "silking", 
                                              "anthesis.silking.interval", 
                                              "grain.number", "grain.yield", "grain.weight"),
                                  selected = "grain.yield"),
                      sliderInput(inputId = "sigma",
                                  label = "Range of interest for sigma:",
                                  min = 2,
                                  max = 25,
                                  value = 1),
                      
                      sliderInput(inputId = "mu",
                                  label = "Range of interest for mu:",
                                  min = -20,
                                  max = 20,
                                  value = 0),
                      
                      sliderInput("enne",
                                  "Number of bins:",
                                  value = 20,
                                  min = 10,
                                  max = 100),
                      plotOutput(outputId = "distPlot")
             ),
             
             # Model summary panel
             tabPanel("Model Summary Viewer", value = "summary",
                      fluidPage(
                        titlePanel("Model"),
                        sidebarLayout(
                          sidebarPanel(
                            selectInput("modelSelect", "Select Model", choices = paste0("mod.", 1:14))
                          ),
                          mainPanel(
                            verbatimTextOutput("modelSummary")
                          )
                        )
                      )
             ),
             
             # Data Visualization panel
             tabPanel("Data Visualization", value = "visualization",
                      fluidPage(
                        titlePanel("Create Data Visualizations"),
                        sidebarLayout(
                          sidebarPanel(
                            selectInput("xvar", "X-axis Variable:", choices = names(df_tesi)),
                            selectInput("yvar", "Y-axis Variable:", choices = names(df_tesi)),
                            selectInput("plotType", "Plot Type:",
                                        choices = c("Scatter Plot" = "scatter", "Histogram Plot" = "histogram" )),
                            actionButton("plotButton", "Generate Plot")
                          ),
                          mainPanel(
                            plotOutput("plotOutput")
                          )
                        )
                      )
             ),
 
             
             # Send Notes panel
             tabPanel("Send Notes", value = "notes",
                      fluidPage(
                        titlePanel("Send Your Notes, Feedback or Suggestions for improvement"),
                        sidebarLayout(
                          sidebarPanel(
                            textAreaInput("notes", "Your Notes:", "", width = '100%', height = '200px'),
                            actionButton("sendEmail", "Send Email"),
                            fileInput("file1", "Choose File to Upload")
                          ),
                          mainPanel(
                            uiOutput("sendEmail"),
                            verbatimTextOutput("upload")
                          )
                        )
                      )
             ),
             
             #menu
             navbarMenu("Menu...",
                        tabPanel("About", value = "about",
                                 p("This application was made with the use of", strong("Shiny!"), p("This study specifically focuses on the grain yield node within a pre-established Bayesian network (see Valleggi et al., 2024) in the domain of agronomical studies. Drawing from Millet et al.'s (2019) genome-wide association study encompassing 256 Maize Zea L. varieties within the DROPS project ((DROught-tolerant yielding PlantS, one of the European Union projects representing a collaborative endeavour poised to advance our comprehension of drought tolerance in crops and provide practical solutions for enhancing yield resilience in challenging environmental conditions), this research aims to investigate relationships in crop production. In fact, leveraging Millet et al.'s groundwork, a rigorous data filtering technique is implemented to ensure dataset completeness, with a final sample size of 143 observations and 21 variables collected from three specific sites (Gaillac, Debrecen and Karlsruhe). Beyond mere replication, this study exclusively centers on the variable grain yield, Bayesian linear models and prior distributions. The research highlights the pivotal role of parameter learning in constructing effective models, addressing the complex probabilistic dependencies inherent in crop production. After scrutinizing 14 models and facing uncertainty between mod.14 and mod.3, the finally selected model is mod.3, deemed the most suitable for performing a detailed characterization. Through statistical methodologies, the study aspires to provide nuanced insights into the relationships between genetic factors, plant performance, and environmental influences. The final aim is to contribute to the knowledge base in crop science by developing a class of Bayesian models able to jointly consider several varieties at different sites, and by exploring the role of prior distributions required in the computation. The analysis of the residuals found no evidence of violations in the assumptions. Model selection has been performed with state-of-the-art leave-one-out algorithms. This structured abstract encapsulates the research's objectives, methodology, and anticipated contributions, providing a comprehensive overview of the study's significance and potential impact in the field."),
                                   tags$img(height = 1000, width = 800, src = "maize.jpg") 
                                 )),
                        
                        tabPanel("Reference", value = "reference",
                                 p("For more details see:", br(),
                                   strong("thesis:'Bayesian Model Selection for Heterogeneous Agronomic Data'. Diabou Mbaye(2024)."),
                                   br(), br(),
                                   strong("BIBLIOGRAPHY"), br(),
                                   tags$ol(
                                     tags$li("Blackwell, David. ‘Basic Statistics’. McGraw Hill (1969)."),
                                     tags$li("Bürkner, Paul-Christian. ‘Brms : An R Package for Bayesian Multilevel Models Using Stan’. In Journal of Statistical Software vol. 80(1). http://www.jstatsoft.org/v80/i01/ (2017)."),
                                     tags$li("Chandrasekaran, B, Annadurai K. and Somasundaram E., ‘A Textbook of Agronomy’. In New Age International (P) Ltd., Publishers (2010)."),
                                     tags$li("Consonni, Guido, Dimitris Fouskakis, Brunero Liseo and Ioannis Ntzoufras. ‘Prior Distributions for Objective Bayesian Analysis’. In Bayesian Analysis vol. 13(2). https://projecteuclid.org/journals/bayesian-analysis/volume-13/issue-2/Prior-Distributions-for-Objective-Bayesian-Analysis/10.1214/18-BA1103.full (2018)."),
                                     tags$li("Debnath, S. C. and Sarkar K. R. ‘Combining ability analysis of grain yield and some of its attributes in maize (zea mays L.)’. Division of Genetics, Indian Agricultural Research Institute (1989)."),
                                     tags$li("‘Food and Agriculture: Key to Achieving the 2030 Agenda for Sustainable Development’. Food and agriculture (FAO). (2016)."),
                                     tags$li("Gelman, Andrew and Jennifer Hill. ‘Data Analysis Using Regression and Multilevel/Hierarchical Models’. Cambridge University Press. (2006)."),
                                     tags$li("Gelman, Andrew, Carlin John B., Stern Hal S., Dunson David B., Vehtari Aki and Rubin Donald B. ‘Bayesian Data Analysis Third Edition (with Errors Fixed as of 15 February 2021)’. (2021)."),
                                     tags$li("Gelman, Andrew, Brooks Steve, Jones L. Galin and Meng Xiao-Li. ‘Handbook of Markov Chain Monte Carlo’. Chapman & Hall. (CRC). (2011)."),
                                     tags$li("Green, K. Highly Structured Stochastic Systems. Canada: John Wiley & Sons Canada, Limited. (2001)."),
                                     tags$li("Kutner, M. H., Nachtsheim, C. J., Neter, J., and Li, W. Applied Linear Statistical Models (5th ed.). McGraw-Hill Irwin. (2004)."),
                                     tags$li("Larson, Martin G. ‘Descriptive Statistics and Graphical Displays’. In Circulation vol. 114(1): 76–81. https://www.ahajournals.org/doi/full/10.1161/CIRCULATIONAHA.105.584474 (2006)."),
                                     tags$li("Mercer, KL, Martınez-Vasquez A and Perales HR. ‘Asymmetrical local adaptation of maize landraces along an altitudinal gradient’. In Evol. Appl. Vol. 1:489–500 (2008)."),
                                     tags$li("Mercer, KL and Perales HR. ‘Structure of local adaptation across the landscape: flowering time and fitness in Mexican maize (Zea mays L. subsp. mays) landraces’. In Genet. Resour. Crop Evol. Vol. 66:27–45 (2019)."),
                                     tags$li("Millet, E. J., Pommier C., Buy M., Nagel A., Kruijer W., Welz-Bolduan T., Lopez J., Richard C., Racz F., Tanzi F., Spitkot T., Canè M.A., Negro S.S., Coupel-Ledru A., Nicolas S.D., Palaffre C., Bauland C., Praud S., Ranc N., Presterl T., Bedo Z., Tuberosa R., Usadel B., Charcosset A., van Eeuwijk F. A., Draye X., Tardieu F. and Welcker C. ‘A Multi-Site Experiment in a Network of European Fields for Assessing the Maize Yield Response to Environmental Scenarios’. doi:10.15454/IASSTN (2019)."),
                                     tags$li("Pernkopf, Franz, Robert Peharz, and Sebastian Tschiatschek. ‘Chapter 18 - Introduction to Probabilistic Graphical Models’. In Academic Press Library in Signal Processing, Academic Press Library in Signal Processing: Vol. 1, eds. Paulo S. R. Diniz, Johan A. K. Suykens, Rama Chellappa, and Sergios Theodoridis. Elsevier, 989–1064. https://www.sciencedirect.com/science/article/pii/B9780123965028000188 (2014)."),
                                     tags$li("Rogers, Anna R. et al. ‘The Importance of Dominance and Genotype-by-Environment Interactions on Grain Yield Variation in a Large-Scale Public Cooperative Maize Experiment’. Oxford. https://academic.oup.com/g3journal/article/doi/10.1093/g3journal/jkaa050/6062399 (2021)."),
                                     tags$li("Romay, MC, Millard MJ, Glaubitz JC, Peiffer JA, Swarts KL, et al. ‘Comprehensive genotyping of the USA national maize inbred seed bank’. In Genome Biol. Vol. 14: R55 (2013)."),
                                     tags$li("Romero, Navarro JA, Willcox M, Burgueno J, Romay MC, Swarts KL, et al. ‘A study of allelic diversity underlying flowering-time adaptation in maize landraces’. In Nat. Genet. Vol. 49:476–480 (2017)."),
                                     tags$li("Ruiz, Corral JA, Dura´ n Puga N, Sanchez Gonzalez JDJ, Ron Parra J, Gonzalez Eguiarte DR, et al. ‘Climatic adaptation and ecological descriptors of 42 Mexican maize races’. In Crop Sci. vol. 48: 1502–1512 (2008)."),
                                     tags$li("Su, Xiaogang, Xin Yan, and Chih-Ling Tsai. ‘Linear Regression’. In WIREs Computational Statistics vol. 4 (3): 275–94. https://onlinelibrary.wiley.com/doi/abs/10.1002/wics.1198 (2012)."),
                                     tags$li("Valleggi, Lorenzo, Scutari Marco and Stefanini M. Federico ‘Learning Bayesian Networks with Heterogeneous Agronomic Data Sets via Mixed-Effect Models and Hierarchical Clustering’ (2024). Engineering Applications of Artificial Intelligence."),
                                     tags$li("Vehtari, A., Gelman, A., & Gabry, J. (2017). Practical Bayesian model evaluation using leave-one-out cross-validation and WAIC. Statistics and Computing, 27(5), 1413–1432. https://doi.org/10.1007/s11222-016-9696-4")
                                   ),
                                   br(),
                                   strong("SITOGRAPHY"), br(),
                                   tags$ul(
                                     tags$li(tags$a(href = "https://www.fao.org/about/about-fao/en/", "Food and Agriculture Organization of the United Nations (FAO)")),
                                     tags$li(tags$a(href = "https://www.inrae.fr/en/about-us", "INRAE - National Research Institute for Agriculture, Food and Environment")),
                                     tags$li(tags$a(href = "https://ec.europa.eu/eurostat/statistics-explained/index.php?title=Agricultural_production_-_crops", "Eurostat - Agricultural Production - Crops")),
                                     tags$li(tags$a(href = "https://cordis.europa.eu/project/id/244374", "DROPS Project - European Commission")),
                                     tags$li(tags$a(href = "https://www.r-project.org/", "R Project for Statistical Computing")),
                                     tags$li(tags$a(href = "https://mc-stan.org/docs/stan-users-guide/index.html", "Stan User’s Guide"))
                                   ),
                                   br(),
                                   strong("For additional resources on Shiny:"),
                                   br(),
                                   tags$a(href = "https://shiny.rstudio.com/reference/shiny/latest/navbarPage.html", "Shiny - navbarPage"), br(),
                                   tags$a(href = "https://shiny.posit.co/r/articles/build/layout-guide/", "Shiny Layout Guide"), br(),
                                   tags$a(href = "https://shiny.posit.co/r/getstarted/shiny-basics/lesson1/index.html", "Getting Started with Shiny")
                                 )
                        )
             )
  )
)

server <- function(input, output, session) {
  
  # Create reactive data frame for selected variables
  variables_selected <- reactive({
    df_tesi %>% select(all_of(input$selected_var))
  })
  
  # Create data table
  output$filetable <- renderDataTable({
    req(input$selected_var)
    datatable(
      data = variables_selected(),
      options = list(pageLength = 10),
      rownames = FALSE
    )
  })
  
  output$trialTable <- renderTable({
    trial_data 
  })
  
  output$dataSummary <- renderPrint({
    summary(df_tesi)
  })
  
  # Download file
  output$download_data <- downloadHandler(
    filename = function() {
      paste0("data.", input$filetype)
    },
    content = function(file) {
      selected_data <- variables_selected()
      if (input$filetype == "csv") {
        write_csv(selected_data, file)
      } else if (input$filetype == "tsv") {
        write_tsv(selected_data, file)
      } else if (input$filetype == "txt") {
        write.table(selected_data, file, sep = "\t", row.names = FALSE)
      } else if (input$filetype == "docx") {
        # Use the 'officer' package to write a Word document
        library(officer)
        doc <- read_docx()
        doc <- body_add_table(doc, value = selected_data, style = "table_template")
        print(doc, target = file)
      }
    }
  )
  
  output$table <- renderTable({
    # Filter df_complete based on the number of observations input
    if (input$observations > 100) {
      df_tesi <- df_complete %>%
        filter(Variety_ID %in% c("B73", "F1808", "F7028")) %>%
        filter(Site %in% c("Gaillac", "Karlsruhe", "Debrecen")) %>%
        sample_n(143)  # Sample 143 observations for df_tesi
    } else {
      df_tesi <- df_complete %>%
        filter(Variety_ID %in% c("B73", "F1808")) %>%
        filter(Site %in% c("Gaillac", "Karlsruhe")) %>%
        sample_n(143)  # Sample 143 observations for df_tesi
    }
    return(df_tesi)
  })
  

  
  output$distPlot <- renderPlot({
    # Estraction of selected variable
    centerDatTB <- tibble(pippo = wDF[[input$var]]) %>%
      mutate(pippo = (pippo - mean(pippo)) * 10 / sd(pippo))
    
    # Creation of theoretic grid
    xgrid <- seq(input$mu - input$sigma * 4, input$mu + input$sigma * 4, length.out = 500)
    
    
    ygrid <- switch(input$dist,
                    "norm" = dnorm(xgrid, mean = input$mu, sd = input$sigma),
                    "unif" = dunif(xgrid, min = input$mu, max = input$sigma),
                    "lnorm" = dlnorm(xgrid, meanlog = input$mu, sdlog = input$sigma),
                    "exp" = dexp(xgrid, rate = 1 / input$sigma))
    
    # Plot + fitted line
    ggplot() +
      geom_line(aes(x = xgrid, y = ygrid), colour = "blue") +  # Linea della distribuzione teorica
      geom_histogram(data = centerDatTB, aes(y = after_stat(density), x = pippo),
                     alpha = 0.3, bins = input$enne, colour = "#8B5A00", fill = "#CD8500") +
      xlab(input$var) +
      theme(axis.text = element_text(size = 12),
            axis.title = element_text(size = 14, face = "plain")) +
      xlim(c(-100, 100))
  })
  
  
  
  # Reactive expression to retrieve and display the selected model's summary
  output$modelSummary <- renderPrint({
    req(input$modelSelect)
    
    # Get the selected model dynamically
    selectedModel <- get(input$modelSelect)
    
    # Display summary of the selected model
    summary(selectedModel)
  })

  
  
  # Definisci plotReady come una variabile reattiva
  plotReady <- reactiveVal(FALSE)  # Inizializzata a FALSE
  
  # Observe the button click to trigger the plot generation
  observeEvent(input$plotButton, {
    req(input$xvar, input$yvar, input$plotType)
    
    plotData <- df_tesi %>% select(all_of(c(input$xvar, input$yvar)))
    
    # Check if selected variables are valid
    if (!(input$xvar %in% colnames(df_tesi)) || !(input$yvar %in% colnames(df_tesi))) {
      shiny::showNotification("Selected variable(s) not found in the dataset.", type = "error")
      plotReady(FALSE)  # Reset plot readiness
      return(NULL)  # Stop further processing
    }
    
    # Check if both selected variables are numeric (for scatter plot)
    if (!is.numeric(plotData[[input$xvar]]) && input$plotType == "scatter") {
      shiny::showNotification("Both selected variables must be numeric for a scatter plot. Please choose numeric variables.", type = "error")
      plotReady(FALSE)  # Reset plot readiness
      return(NULL)  # Stop further processing
    }
    
    # Se tutte le condizioni sono soddisfatte, abilita il rendering del grafico
    plotReady(TRUE)
  })
  
  # Render the plot only if plotReady is TRUE
  output$plotOutput <- renderPlot({
    if (!plotReady()) {
      return(NULL)  # Non mostrare nulla se il grafico non è pronto
    }
    
    # Genera il grafico in base al tipo di grafico selezionato
    plotData <- df_tesi %>% select(all_of(c(input$xvar, input$yvar)))
    
    # Scatter plot
    if (input$plotType == "scatter") {
      ggplot(plotData, aes_string(x = input$xvar, y = input$yvar)) + 
        geom_point() +
        labs(title = paste("Scatter plot of", input$xvar, "vs", input$yvar)) +
        theme_minimal()
      
      # Histogram
    } else if (input$plotType == "histogram") {
      if (is.numeric(plotData[[input$xvar]])) {
        ggplot(plotData, aes_string(x = input$xvar)) + 
          geom_histogram(bins = 30) +  # Histogram for numeric variables
          labs(title = paste("Histogram of", input$xvar)) +
          theme_minimal()
      } else {
        ggplot(plotData, aes_string(x = input$xvar)) + 
          geom_bar(stat = "count") +  # Bar plot for discrete variables
          labs(title = paste("Bar plot of", input$xvar)) +
          theme_minimal()
      }
    }
  })
  
  

  
  # Send email
  observeEvent(input$sendEmail, {
    notes <- input$notes
    mailto_link <- sprintf(
      "mailto:diabou22@libero.it?subject=User%%20Notes&body=%s",
      URLencode(notes)
    )
    output$sendEmail <- renderUI({
      tags$a(href = mailto_link, "Send Email", target = "_blank", class = "btn btn-primary")
    })
  })
  
  # Display uploaded file information
  output$upload <- renderPrint({
    if (!is.null(input$file1)) {
      paste("File name:", input$file1$name, "<br>",
            "File type:", input$file1$type, "<br>",
            "File size:", prettyNum(input$file1$size, big.mark = ","), "bytes")
    }
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
