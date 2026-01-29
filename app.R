# Dengue Variant Tracker - Interactive Dashboard
# Shiny application for exploring dengue virus mutation patterns

library(shiny)
library(shinydashboard)
library(ggplot2)
library(dplyr)
library(plotly)
library(DT)
library(Biostrings)

# ===== Load Data =====
load_data <- function() {
  tryCatch({
    # Load processed analysis results
    qc_data <- read.csv("data/processed/qc_summary.csv", stringsAsFactors = FALSE)
    motif_data <- read.csv("data/processed/motif_matches.csv", stringsAsFactors = FALSE)
    motif_summary <- read.csv("data/processed/motif_summary.csv", stringsAsFactors = FALSE)
    
    # Load sequences for custom motif search
    seqs <- readDNAStringSet("data/processed/filtered_sequences.fasta")
    
    list(
      qc = qc_data,
      motifs = motif_data,
      motif_summary = motif_summary,
      sequences = seqs,
      loaded = TRUE
    )
  }, error = function(e) {
    list(
      loaded = FALSE,
      error_message = paste("Error loading data:", e$message, 
                           "\n\nPlease run qc_analysis.R first.")
    )
  })
}

data <- load_data()

# ===== UI Definition =====
ui <- dashboardPage(
  skin = "blue",
  
  # Header
  dashboardHeader(
    title = "🦟 Dengue Variant Tracker",
    titleWidth = 300
  ),
  
  # Sidebar
  dashboardSidebar(
    width = 300,
    sidebarMenu(
      menuItem("Overview", tabName = "overview", icon = icon("dashboard")),
      menuItem("Sequence Quality", tabName = "quality", icon = icon("check-circle")),
      menuItem("Motif Explorer", tabName = "motifs", icon = icon("search")),
      menuItem("Custom Analysis", tabName = "custom", icon = icon("dna")),
      menuItem("Data Table", tabName = "data", icon = icon("table")),
      menuItem("About", tabName = "about", icon = icon("info-circle"))
    ),
    
    hr(),
    
    div(style = "padding: 15px;",
        h5("Dataset Info", style = "color: white; font-weight: bold;"),
        if (data$loaded) {
          tagList(
            p(paste("Sequences:", nrow(data$qc)), style = "color: white; margin: 5px 0;"),
            p(paste("Motifs Found:", nrow(data$motifs)), style = "color: white; margin: 5px 0;"),
            p("Source: NCBI Virus", style = "color: white; margin: 5px 0; font-size: 11px;")
          )
        } else {
          p("Data not loaded", style = "color: #ff6b6b;")
        }
    )
  ),
  
  # Body
  dashboardBody(
    tags$head(
      tags$style(HTML("
        .content-wrapper { background-color: #f4f6f9; }
        .box { border-radius: 8px; box-shadow: 0 2px 4px rgba(0,0,0,0.1); }
        .info-box { border-radius: 8px; }
        .small-box { border-radius: 8px; }
        .nav-tabs-custom { border-radius: 8px; }
        h2 { color: #2c3e50; font-weight: bold; }
      "))
    ),
    
    tabItems(
      # ===== OVERVIEW TAB =====
      tabItem(
        tabName = "overview",
        if (!data$loaded) {
          box(
            width = 12,
            status = "danger",
            title = "Error Loading Data",
            solidHeader = TRUE,
            p(data$error_message),
            p("Please ensure you have run the data download and analysis scripts:"),
            tags$ol(
              tags$li(code("./download_data.sh")),
              tags$li(code("Rscript qc_analysis.R"))
            )
          )
        } else {
          tagList(
            fluidRow(
              valueBox(
                value = nrow(data$qc),
                subtitle = "Total Sequences",
                icon = icon("dna"),
                color = "blue",
                width = 3
              ),
              valueBox(
                value = format(sum(data$qc$Length), big.mark = ","),
                subtitle = "Total Nucleotides",
                icon = icon("ruler"),
                color = "green",
                width = 3
              ),
              valueBox(
                value = paste0(round(mean(data$qc$GC_Content) * 100, 1), "%"),
                subtitle = "Mean GC Content",
                icon = icon("chart-pie"),
                color = "yellow",
                width = 3
              ),
              valueBox(
                value = nrow(data$motifs),
                subtitle = "Motif Matches",
                icon = icon("bullseye"),
                color = "red",
                width = 3
              )
            ),
            
            fluidRow(
              box(
                title = "Sequence Length Distribution",
                status = "primary",
                solidHeader = TRUE,
                width = 6,
                plotlyOutput("overview_length_plot", height = "300px")
              ),
              box(
                title = "GC Content Distribution",
                status = "success",
                solidHeader = TRUE,
                width = 6,
                plotlyOutput("overview_gc_plot", height = "300px")
              )
            ),
            
            fluidRow(
              box(
                title = "Top Motifs Found",
                status = "warning",
                solidHeader = TRUE,
                width = 12,
                DTOutput("overview_motif_table")
              )
            )
          )
        }
      ),
      
      # ===== QUALITY TAB =====
      tabItem(
        tabName = "quality",
        if (!data$loaded) {
          box(width = 12, status = "danger", title = "Data Not Loaded", 
              p(data$error_message))
        } else {
          tagList(
            fluidRow(
              box(
                title = "Quality Control Metrics",
                status = "info",
                solidHeader = TRUE,
                width = 12,
                tabsetPanel(
                  tabPanel(
                    "Length Analysis",
                    br(),
                    plotlyOutput("qc_length_detail", height = "400px")
                  ),
                  tabPanel(
                    "Base Composition",
                    br(),
                    plotlyOutput("qc_nucleotide_comp", height = "400px")
                  ),
                  tabPanel(
                    "GC vs Length",
                    br(),
                    plotlyOutput("qc_gc_scatter", height = "400px")
                  )
                )
              )
            ),
            
            fluidRow(
              box(
                title = "Summary Statistics",
                status = "primary",
                solidHeader = TRUE,
                width = 12,
                DTOutput("qc_stats_table")
              )
            )
          )
        }
      ),
      
      # ===== MOTIF EXPLORER TAB =====
      tabItem(
        tabName = "motifs",
        if (!data$loaded) {
          box(width = 12, status = "danger", title = "Data Not Loaded", 
              p(data$error_message))
        } else {
          tagList(
            fluidRow(
              box(
                title = "Filter Motifs",
                status = "primary",
                solidHeader = TRUE,
                width = 3,
                selectInput(
                  "motif_filter",
                  "Select Motif:",
                  choices = c("All", unique(data$motifs$Motif)),
                  selected = "All"
                ),
                hr(),
                checkboxInput("show_positions", "Show Position Details", value = FALSE),
                actionButton("refresh_motif", "Refresh", icon = icon("sync"), 
                            class = "btn-primary btn-block")
              ),
              
              box(
                title = "Motif Frequency",
                status = "success",
                solidHeader = TRUE,
                width = 9,
                plotlyOutput("motif_frequency_plot", height = "400px")
              )
            ),
            
            fluidRow(
              box(
                title = "Motif Match Details",
                status = "warning",
                solidHeader = TRUE,
                width = 12,
                DTOutput("motif_detail_table")
              )
            )
          )
        }
      ),
      
      # ===== CUSTOM ANALYSIS TAB =====
      tabItem(
        tabName = "custom",
        if (!data$loaded) {
          box(width = 12, status = "danger", title = "Data Not Loaded", 
              p(data$error_message))
        } else {
          tagList(
            fluidRow(
              box(
                title = "Custom Motif Search",
                status = "primary",
                solidHeader = TRUE,
                width = 12,
                p("Search for custom DNA motifs in the dengue sequences. Enter a DNA sequence pattern using standard nucleotide codes."),
                
                fluidRow(
                  column(
                    6,
                    textInput(
                      "custom_motif",
                      "DNA Motif Pattern:",
                      value = "ATG",
                      placeholder = "e.g., ATG, CACAG, AATAAA"
                    ),
                    helpText("Use A, T, G, C, or IUPAC ambiguity codes (N, R, Y, etc.)")
                  ),
                  column(
                    6,
                    numericInput(
                      "max_mismatch",
                      "Maximum Mismatches:",
                      value = 0,
                      min = 0,
                      max = 3
                    ),
                    helpText("Allow up to N mismatches for fuzzy matching")
                  )
                ),
                
                actionButton("search_custom", "Search Motif", 
                            icon = icon("search"), 
                            class = "btn-success btn-lg")
              )
            ),
            
            fluidRow(
              box(
                title = "Search Results",
                status = "success",
                solidHeader = TRUE,
                width = 12,
                uiOutput("custom_results_summary"),
                hr(),
                plotlyOutput("custom_motif_plot", height = "350px"),
                hr(),
                DTOutput("custom_results_table")
              )
            )
          )
        }
      ),
      
      # ===== DATA TABLE TAB =====
      tabItem(
        tabName = "data",
        if (!data$loaded) {
          box(width = 12, status = "danger", title = "Data Not Loaded", 
              p(data$error_message))
        } else {
          box(
            title = "Complete Dataset",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            DTOutput("full_data_table")
          )
        }
      ),
      
      # ===== ABOUT TAB =====
      tabItem(
        tabName = "about",
        fluidRow(
          box(
            title = "About This Dashboard",
            status = "info",
            solidHeader = TRUE,
            width = 12,
            h3("Dengue Variant Tracker"),
            p("An interactive bioinformatics dashboard for exploring dengue virus genomic variants 
              and mutation patterns using R and Bioconductor."),
            
            h4("Features:"),
            tags$ul(
              tags$li("Quality control analysis of dengue virus sequences"),
              tags$li("Automated motif detection for known mutation sites"),
              tags$li("Interactive visualizations of sequence characteristics"),
              tags$li("Custom motif search capabilities"),
              tags$li("Export-ready data tables")
            ),
            
            h4("Data Source:"),
            p("Sequences obtained from the NCBI Virus Database (public domain)."),
            tags$a(href = "https://www.ncbi.nlm.nih.gov/labs/virus/", 
                  "Visit NCBI Virus Database", target = "_blank"),
            
            h4("Technologies:"),
            p("R, Bioconductor (Biostrings, ShortRead), Shiny, ggplot2, plotly, dplyr"),
            
            h4("Ethical Considerations:"),
            tags$ul(
              tags$li("Uses only public, anonymized viral sequence data"),
              tags$li("No personal health information is collected or displayed"),
              tags$li("All data sources are properly cited"),
              tags$li("Code is open-source and reproducible")
            ),
            
            hr(),
            p("Developed for bioinformatics portfolio demonstration.", 
              style = "font-style: italic; color: #7f8c8d;"),
            p("For questions or feedback, please refer to the GitHub repository.", 
              style = "font-style: italic; color: #7f8c8d;")
          )
        )
      )
    )
  )
)

# ===== Server Logic =====
server <- function(input, output, session) {
  
  # Overview tab plots
  output$overview_length_plot <- renderPlotly({
    if (!data$loaded) return(NULL)
    
    p <- ggplot(data$qc, aes(x = Length)) +
      geom_histogram(bins = 30, fill = "#3498db", alpha = 0.7) +
      labs(x = "Sequence Length (bp)", y = "Count") +
      theme_minimal()
    
    ggplotly(p) %>% layout(hovermode = "x unified")
  })
  
  output$overview_gc_plot <- renderPlotly({
    if (!data$loaded) return(NULL)
    
    p <- ggplot(data$qc, aes(x = GC_Content * 100)) +
      geom_histogram(bins = 20, fill = "#2ecc71", alpha = 0.7) +
      labs(x = "GC Content (%)", y = "Count") +
      theme_minimal()
    
    ggplotly(p) %>% layout(hovermode = "x unified")
  })
  
  output$overview_motif_table <- renderDT({
    if (!data$loaded) return(NULL)
    
    datatable(
      data$motif_summary %>% select(Motif, Description, Total_Occurrences, 
                                    Sequences_With_Motif),
      options = list(pageLength = 5, dom = 't'),
      rownames = FALSE,
      colnames = c("Motif", "Description", "Total Occurrences", "Sequences")
    )
  })
  
  # Quality tab plots
  output$qc_length_detail <- renderPlotly({
    if (!data$loaded) return(NULL)
    
    p <- ggplot(data$qc, aes(x = reorder(Sequence_ID, Length), y = Length)) +
      geom_col(fill = "#9b59b6", alpha = 0.8) +
      labs(x = "Sequence", y = "Length (bp)", title = "Individual Sequence Lengths") +
      theme_minimal() +
      theme(axis.text.x = element_blank())
    
    ggplotly(p)
  })
  
  output$qc_nucleotide_comp <- renderPlotly({
    if (!data$loaded) return(NULL)
    
    nuc_data <- data.frame(
      Nucleotide = rep(c("A", "T", "G", "C"), each = nrow(data$qc)),
      Count = c(data$qc$A_Count, data$qc$T_Count, data$qc$G_Count, data$qc$C_Count),
      Sequence = rep(data$qc$Sequence_ID, 4)
    )
    
    p <- ggplot(nuc_data, aes(x = Nucleotide, y = Count, fill = Nucleotide)) +
      geom_boxplot(alpha = 0.7) +
      scale_fill_manual(values = c(A = "#e74c3c", T = "#3498db", 
                                   G = "#f39c12", C = "#1abc9c")) +
      labs(title = "Nucleotide Count Distribution", y = "Count") +
      theme_minimal() +
      theme(legend.position = "none")
    
    ggplotly(p)
  })
  
  output$qc_gc_scatter <- renderPlotly({
    if (!data$loaded) return(NULL)
    
    p <- ggplot(data$qc, aes(x = Length, y = GC_Content * 100, 
                             text = Sequence_ID)) +
      geom_point(color = "#e67e22", size = 3, alpha = 0.6) +
      geom_smooth(method = "lm", se = TRUE, color = "#34495e") +
      labs(x = "Sequence Length (bp)", y = "GC Content (%)", 
           title = "GC Content vs Sequence Length") +
      theme_minimal()
    
    ggplotly(p, tooltip = c("text", "x", "y"))
  })
  
  output$qc_stats_table <- renderDT({
    if (!data$loaded) return(NULL)
    
    summary_stats <- data.frame(
      Metric = c("Mean Length", "Median Length", "Min Length", "Max Length",
                "Mean GC%", "Median GC%", "Min GC%", "Max GC%"),
      Value = c(
        round(mean(data$qc$Length), 0),
        round(median(data$qc$Length), 0),
        min(data$qc$Length),
        max(data$qc$Length),
        round(mean(data$qc$GC_Content) * 100, 2),
        round(median(data$qc$GC_Content) * 100, 2),
        round(min(data$qc$GC_Content) * 100, 2),
        round(max(data$qc$GC_Content) * 100, 2)
      )
    )
    
    datatable(summary_stats, options = list(dom = 't'), rownames = FALSE)
  })
  
  # Motif explorer
  filtered_motifs <- reactive({
    if (!data$loaded) return(NULL)
    
    input$refresh_motif  # Dependency for refresh button
    
    if (input$motif_filter == "All") {
      data$motifs
    } else {
      data$motifs %>% filter(Motif == input$motif_filter)
    }
  })
  
  output$motif_frequency_plot <- renderPlotly({
    if (!data$loaded) return(NULL)
    
    motif_counts <- filtered_motifs() %>%
      group_by(Motif) %>%
      summarise(Count = n(), .groups = 'drop')
    
    p <- ggplot(motif_counts, aes(x = reorder(Motif, Count), y = Count)) +
      geom_col(fill = "#16a085", alpha = 0.8) +
      coord_flip() +
      labs(x = "Motif", y = "Occurrences", title = "Motif Match Frequency") +
      theme_minimal()
    
    ggplotly(p)
  })
  
  output$motif_detail_table <- renderDT({
    if (!data$loaded) return(NULL)
    
    display_data <- if (input$show_positions) {
      filtered_motifs()
    } else {
      filtered_motifs() %>% 
        select(Sequence_ID, Motif, Description, Count) %>%
        distinct()
    }
    
    datatable(
      display_data,
      options = list(pageLength = 10, scrollX = TRUE),
      rownames = FALSE,
      filter = 'top'
    )
  })
  
  # Custom motif search
  custom_search_results <- eventReactive(input$search_custom, {
    if (!data$loaded) return(NULL)
    
    motif_pattern <- toupper(input$custom_motif)
    max_mm <- input$max_mismatch
    
    withProgress(message = 'Searching sequences...', value = 0, {
      tryCatch({
        # Search for motif
        matches <- vmatchPattern(
          motif_pattern, 
          data$sequences,
          max.mismatch = max_mm
        )
        
        # Process results
        results_list <- lapply(seq_along(matches), function(i) {
          incProgress(1 / length(matches))
          if (length(matches[[i]]) > 0) {
            data.frame(
              Sequence_ID = names(data$sequences)[i],
              Motif = motif_pattern,
              Position = start(matches[[i]]),
              Match_Count = length(matches[[i]]),
              stringsAsFactors = FALSE
            )
          } else {
            NULL
          }
        })
        
        results <- do.call(rbind, results_list)
        
        if (is.null(results) || nrow(results) == 0) {
          list(
            found = FALSE,
            message = paste("No matches found for motif:", motif_pattern)
          )
        } else {
          list(
            found = TRUE,
            data = results,
            total_matches = nrow(results),
            sequences_matched = n_distinct(results$Sequence_ID)
          )
        }
      }, error = function(e) {
        list(
          found = FALSE,
          message = paste("Error during search:", e$message)
        )
      })
    })
  })
  
  output$custom_results_summary <- renderUI({
    results <- custom_search_results()
    if (is.null(results)) return(NULL)
    
    if (!results$found) {
      div(
        class = "alert alert-warning",
        icon("exclamation-triangle"),
        " ",
        results$message
      )
    } else {
      div(
        class = "alert alert-success",
        icon("check-circle"),
        sprintf(" Found %d matches in %d sequences", 
                results$total_matches, results$sequences_matched)
      )
    }
  })
  
  output$custom_motif_plot <- renderPlotly({
    results <- custom_search_results()
    if (is.null(results) || !results$found) return(NULL)
    
    match_summary <- results$data %>%
      group_by(Sequence_ID) %>%
      summarise(Matches = n(), .groups = 'drop')
    
    p <- ggplot(match_summary, aes(x = reorder(Sequence_ID, Matches), y = Matches)) +
      geom_col(fill = "#c0392b", alpha = 0.8) +
      coord_flip() +
      labs(x = "Sequence", y = "Number of Matches", 
           title = paste("Matches for", input$custom_motif)) +
      theme_minimal() +
      theme(axis.text.y = element_text(size = 8))
    
    ggplotly(p)
  })
  
  output$custom_results_table <- renderDT({
    results <- custom_search_results()
    if (is.null(results) || !results$found) return(NULL)
    
    datatable(
      results$data,
      options = list(pageLength = 10, scrollX = TRUE),
      rownames = FALSE,
      filter = 'top'
    )
  })
  
  # Data table tab
  output$full_data_table <- renderDT({
    if (!data$loaded) return(NULL)
    
    datatable(
      data$qc,
      options = list(pageLength = 25, scrollX = TRUE),
      rownames = FALSE,
      filter = 'top'
    )
  })
}

# Run the application
shinyApp(ui = ui, server = server)
