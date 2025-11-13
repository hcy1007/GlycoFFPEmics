library(shiny)
library(shinydashboard)
library(dplyr)
library(tidyr)
library(ggplot2)
library(VIM)
library(ggcorrplot)
library(ggfortify)
library(ggalluvial)
library(DT)
library(shinyjs)
library(stringr)
library(purrr)
library(pheatmap)
library(plotly)
library(htmlwidgets)
library(corrplot)
library(impute)
library(ggvenn)
library(ggVennDiagram)
library(factoextra)
library(MASS)
library(e1071)
library(clusterProfiler)
library(org.Hs.eg.db)
library(patchwork)
library(ggrepel)
library(ggpubr)
library(base64enc)

# Set maximum upload size to 1000MB
options(shiny.maxRequestSize = 1000*1024^2)

# Load helper functions
source("helper_functions.R")

# Add resource path for images
addResourcePath("images", "www")
addResourcePath("www", "www")
addResourcePath("assets", "assets")

# Define UI
ui <- dashboardPage(
  dashboardHeader(
    title = tags$div(
      style = "display: flex; align-items: center; justify-content: center; height: 50px; width: 200px; font-size: 18px; font-weight: bold; color: #134E4A;",
      "GlycoFFPEmics"
    ),
    titleWidth = 300
  ),
  
  dashboardSidebar(
    sidebarMenu(
      id = "tabs",
      menuItem("Home", tabName = "main_page", icon = icon("home")),
      menuItem("Data Processing", tabName = "data_processing", icon = icon("database")),
      menuItem("Data Analysis", tabName = "data_analysis", icon = icon("chart-line"),
        menuSubItem("Data Upload", tabName = "data_upload", icon = icon("upload")),
        menuSubItem("Data Preview and Quality Control", tabName = "data_preview", icon = icon("table")),
        menuSubItem("Glycomics Analysis", tabName = "glycomics", icon = icon("dna")),
        menuSubItem("Glycoproteomics Qualitative Analysis", tabName = "glycopeptide", icon = icon("microscope")),
        menuSubItem("Glycoproteomics Quantitative Analysis", tabName = "igp_quant", icon = icon("chart-bar")),
        menuSubItem("Export Results", tabName = "export", icon = icon("download")),
        menuSubItem("Help", tabName = "help", icon = icon("question-circle"))
      )
    ),
    width = 250
  ),
  
  dashboardBody(
    tags$head(
      tags$link(rel = "stylesheet", type = "text/css", href = "custom_sidebar.css"),
      tags$link(rel = "stylesheet", type = "text/css", href = "custom_theme.css"),
      tags$style(HTML("
        .main-page-btn {
          transition: all 0.3s ease;
          border-radius: 8px;
          font-weight: bold;
          text-transform: uppercase;
          letter-spacing: 1px;
        }
        .main-page-btn:hover {
          transform: translateY(-2px);
          box-shadow: 0 4px 8px rgba(0,0,0,0.2);
        }
        .main-page-box {
          transition: all 0.3s ease;
          border-radius: 10px;
        }
        .main-page-box:hover {
          transform: translateY(-5px);
          box-shadow: 0 8px 16px rgba(0,0,0,0.1);
        }
        .welcome-title {
          background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
          -webkit-background-clip: text;
          -webkit-text-fill-color: transparent;
          background-clip: text;
          font-weight: bold;
        }
      "))
    ),
    useShinyjs(),
    tabItems(
      # Home page
      tabItem(tabName = "main_page",
              # Large empty space to show background logo
              div(style = "height: 80vh;"),
              fluidRow(
                column(width = 6,
                       box(width = 12, 
                           title = "Data Processing", 
                           status = "primary",
                           solidHeader = TRUE,
                           height = 300,
                           class = "main-page-box",
                           div(style = "text-align: center; padding: 20px;",
                               icon("database", size = "4x", style = "color: #3498db; margin-bottom: 20px;"),
                               h3("Data Processing Module"),
                               p("For data preprocessing, cleaning, normalization and other operations"),
                               br(),
                               actionButton("go_data_processing", "Enter Data Processing", 
                                          class = "btn btn-primary btn-lg main-page-btn",
                                          style = "width: 200px;")
                           )
                       )
                ),
                column(width = 6,
                       box(width = 12, 
                           title = "Data Analysis", 
                           status = "success",
                           solidHeader = TRUE,
                           height = 300,
                           class = "main-page-box",
                           div(style = "text-align: center; padding: 20px;",
                               icon("chart-line", size = "4x", style = "color: #27ae60; margin-bottom: 20px;"),
                               h3("Data Analysis Module"),
                               p("Includes glycomics analysis, glycopeptide qualitative analysis, glycopeptide quantitative analysis and other functions"),
                               br(),
                               actionButton("go_data_analysis", "Enter Data Analysis", 
                                          class = "btn btn-success btn-lg main-page-btn",
                                          style = "width: 200px;")
                           )
                       )
                )
              ),
              fluidRow(
                column(width = 12,
                       box(width = 12,
                           title = "Platform Introduction",
                           status = "info",
                           solidHeader = TRUE,
                           p("GlycoFFPEmics is a comprehensive platform specifically designed for glycomics and glycoproteomics data analysis."),
                           p("Main functions include:"),
                           tags$ul(
                             tags$li("Data Processing: Data preprocessing, cleaning, normalization"),
                             tags$li("Glycomics Analysis: Glycan composition analysis, differential expression analysis"),
                             tags$li("Glycopeptide Qualitative Analysis: Glycopeptide identification, structural analysis"),
                             tags$li("Glycopeptide Quantitative Analysis: Expression analysis, differential analysis, functional annotation")
                           )
                       )
                )
              )
      ),
      
      # Data processing page
      tabItem(tabName = "data_processing",
              fluidRow(
                box(width = 12,
                    title = "Batch Upload Raw CSV Files",
                    fileInput("raw_csv_files", "Select Multiple Raw Data Files (csv)", multiple = TRUE, accept = c(".csv")),
                    actionButton("run_processing", "Start Data Processing", icon = icon("play")),
                    helpText("Please upload raw glycopeptide quantitative csv files. Click the button to automatically complete intra-group and inter-group normalization.")
                )
              ),
              fluidRow(
                box(width = 12,
                    title = "Normalization Results Preview",
                    DTOutput("processing_result_table"),
                    downloadButton("download_processing_result", "Download Normalization Results")
                )
              ),
              fluidRow(
                box(width = 6,
                    title = "Normalized Intensity Distribution Histogram",
                    plotOutput("processing_hist")
                ),
                box(width = 6,
                    title = "Normalized Intensity Boxplot",
                    plotOutput("processing_boxplot")
                )
              )
      ),
      
      # Data upload page
      tabItem(tabName = "data_upload",
              fluidRow(
                box(width = 12,
                    title = "Data File Upload",
                    fileInput("zujian_data", "Select Group Normalization Data File (.csv)",
                              accept = c("text/csv", ".csv")),
                    fileInput("patient_info", "Select Patient Information File (.csv)",
                              accept = c("text/csv", ".csv")),
                    fileInput("glycan_type", "Select Glycan Type File (.csv)",
                              accept = c("text/csv", ".csv")),
                    helpText("Please ensure the uploaded file format is correct, refer to the \"Help\" page for details.")
                )
              ),
              fluidRow(
                box(width = 12,
                    title = "Data Upload Status",
                    verbatimTextOutput("upload_status")
                )
              )
      ),
      
      # Data preview and quality control page
      tabItem(tabName = "data_preview",
              fluidRow(
                box(width = 12,
                    title = "Data Preview",
                    tabsetPanel(
                      tabPanel("Group Normalization Data", DTOutput("zujian_data_preview")),
                      tabPanel("Patient Information Data", DTOutput("patient_info_preview")),
                      tabPanel("Glycan Type Data", DTOutput("glycan_type_preview"))
                    )
                )
              ),
              fluidRow(
                box(width = 12,
                    title = "Data Visualization and Quality Control",
                    tabsetPanel(
                      tabPanel("Sample Matching Status", verbatimTextOutput("sample_matching_status")),
                      tabPanel("Intensity Distribution", plotOutput("intensity_distribution_plot")),
                      tabPanel("Missing Value Overview", plotOutput("missing_value_plot"))
                    )
                )
              )
      ),
      
      # Glycomics analysis page
      tabItem(tabName = "glycomics",
              fluidRow(
                box(width = 12,
                    title = "Glycomics Analysis Parameters",
                    numericInput("valid_ratio", "Effective Value Proportion Threshold:", 0.5, min = 0, max = 1),
                    actionButton("run_glycomics", "Run Glycomics Analysis")
                ),
                box(width = 12,
                    title = "Glycomics Analysis Results",
                    tabsetPanel(
                      tabPanel("Molar Mass Distribution",
                               plotOutput("glycomics_mass_plot"),
                               downloadButton("download_mass_plot", "Download Image")),
                      tabPanel("Molar Mass Distribution (with Significant Annotation)",
                               plotOutput("glycomics_mass_sig_plot"),
                               downloadButton("download_mass_sig_plot", "Download Image")),
                      tabPanel("Significant Difference Volcano Plot",
                               plotOutput("glycomics_significant_plot"),
                               downloadButton("download_significant_plot", "Download Image")),
                      tabPanel("t-test Results",
                               DTOutput("glycomics_ttest_table"),
                               downloadButton("download_glycomics_ttest", "Download t-test Results")),
                      tabPanel("Normalized Data",
                               DTOutput("glycomics_table"),
                               downloadButton("download_glycomics_data", "Download Data"))
                    )
                )
              )
      ),
      
      # Glycopeptide qualitative analysis page
      tabItem(tabName = "glycopeptide",
              fluidRow(
                box(width = 12,
                    title = "Glycopeptide Analysis Parameters",
                    actionButton("run_glycopeptide", "Run Glycopeptide Analysis")
                ),
                box(width = 12,
                    title = "Glycopeptide Analysis Results",
                    tabsetPanel(
                      tabPanel("Sample Distribution Boxplot",
                               plotOutput("glycopeptide_boxplot"),
                               downloadButton("download_boxplot", "Download Image")),
                      tabPanel("Correlation Heatmap",
                               plotOutput("correlation_heatmap"),
                               downloadButton("download_heatmap", "Download Image")),
                      tabPanel("PCA Analysis Plot",
                               plotOutput("pca_plot"),
                               downloadButton("download_pca", "Download Image")),
                      tabPanel("Glycosylation Site Distribution",
                               plotOutput("glycosite_distribution"),
                               downloadButton("download_glycosite", "Download Image")),
                      tabPanel("Sankey Diagram",
                               plotOutput("sankey_plot"),
                               downloadButton("download_sankey", "Download Image")),
                      tabPanel("Unique Counts", 
                               plotOutput("unique_counts_plot"),
                               downloadButton("download_unique_counts_plot", "Download Image")),
                      tabPanel("Data Table",
                               DTOutput("glycopeptide_table"),
                               downloadButton("download_glycopeptide_data", "Download Data")),
                      tabPanel("50% Effective Values Retained _ knn Imputation Before",
                               DTOutput("knn_input_table"),
                               downloadButton("download_knn_input", "Download csv")),
                      tabPanel("Sample Repeats Correlation",
                               DTOutput("correlation_matrix_table"),
                               downloadButton("download_correlation_matrix", "Download csv")),
                      tabPanel("Imputation of Missing Values _ zscore After",
                               DTOutput("zscore_table"),
                               downloadButton("download_zscore", "Download csv"))
                    )
                )
              )
      ),
      
      # Glycopeptide quantitative analysis page
      tabItem(tabName = "igp_quant",
        tabBox(width = 12,
          tabPanel("Data Upload",
            fluidRow(
              box(width = 12, title = "Data Upload",
                fileInput("igp_expr_file", "Upload Imputed Missing Values _ zscore After.csv", accept = c(".csv")),
                fileInput("igp_group_file", "Upload sample_group_PCA.csv", accept = c(".csv")),
                helpText("Please upload the quantitative expression matrix and group information file.")
              )
            ),
            fluidRow(
              box(width = 6, title = "Expression Data Preview", DTOutput("igp_expr_preview")),
              box(width = 6, title = "Group Information Preview", DTOutput("igp_group_preview"))
            )
          ),
          tabPanel("Heatmap Analysis",
            fluidRow(
              box(width = 12, title = "Heatmap Analysis",
                actionButton("run_igp_heatmap", "Run Heatmap Analysis")
              )
            ),
            tabBox(width = 12,
              tabPanel("CRLM vs CRC (Clustering)",
                plotOutput("igp_heatmap_crc"),
                downloadButton("download_igp_heatmap_crc", "Download CRLM vs CRC Heatmap")
              ),
              tabPanel("CRLM vs LIHC (Clustering)",
                plotOutput("igp_heatmap_lihc"),
                downloadButton("download_igp_heatmap_lihc", "Download CRLM vs LIHC Heatmap")
              ),
              tabPanel("Merged Gene Selection (Clustering)",
                plotOutput("igp_heatmap_combined"),
                downloadButton("download_igp_heatmap_combined", "Download Combined Heatmap (Clustering)")
              ),
              tabPanel("Merged Gene Selection (Non-clustering Columns)",
                plotOutput("igp_heatmap_combined_nocluster"),
                downloadButton("download_igp_heatmap_combined_nocluster", "Download Combined Heatmap (Non-clustering Columns)")
              ),
              tabPanel("zscore Standardization (Non-clustering Columns)",
                plotOutput("igp_heatmap_zscore_nocluster"),
                downloadButton("download_igp_heatmap_zscore_nocluster", "Download zscore Heatmap (Non-clustering Columns)")
              ),
              tabPanel("zscore Standardization (Clustering Columns)",
                plotOutput("igp_heatmap_zscore_cluster"),
                downloadButton("download_igp_heatmap_zscore_cluster", "Download zscore Heatmap (Clustering Columns)")
              ),
              tabPanel("zscore Standardization (Clustering Columns + Non-clustering Columns)",
                plotOutput("igp_heatmap_zscore"),
                downloadButton("download_igp_heatmap_zscore", "Download zscore Heatmap (Clustering Columns + Non-clustering Columns)")
              )
            )
          ),
          tabPanel("Regional Heatmap Analysis",
            fluidRow(
              box(width = 12, title = "Regional Analysis Input",
                actionButton("run_regional_analysis", "Run Regional Analysis"),
                helpText('This analysis will use the glycan-glycantype.csv file uploaded on the "Data Upload" page, as well as the two files in the "Data Upload" tab.')
              )
            ),
            fluidRow(
              box(width = 12, title = "Analysis Results",
                tabsetPanel(
                  tabPanel("Mean Ratio Scatter Plot",
                    plotOutput("regional_scatter_plot"),
                    downloadButton("download_regional_scatter", "Download Scatter Plot")
                  ),
                  tabPanel("Regional Statistics",
                    DTOutput("regional_counts_table")
                  ),
                  tabPanel("Heatmap of Differences (Non-central Region)",
                    plotOutput("regional_heatmap_diff"),
                    downloadButton("download_regional_heatmap_diff", "Download Heatmap of Differences")
                  ),
                  tabPanel("z-score Heatmap (Non-central Region)",
                    plotOutput("regional_heatmap_zscore"),
                    downloadButton("download_regional_heatmap_zscore", "Download z-score Heatmap")
                  )
                )
              )
            )
          ),
          tabPanel("Venn Diagram Analysis",
            fluidRow(
              box(width = 12, title = "Venn Diagram Analysis Input",
                actionButton("run_venn_analysis", "Run Venn Diagram Analysis"),
                helpText('This analysis will use the glycan-glycantype.csv file uploaded on the "Data Upload" page, as well as the two files in the "Data Upload" tab.')
              )
            ),
            fluidRow(
              box(width = 12, title = "Venn Diagram Analysis Results",
                tabsetPanel(
                  tabPanel("Glycopeptides Determined by Two Methods",
                    plotOutput("venn_plot_1"),
                    downloadButton("download_venn_plot_1", "Download Venn Diagram 1")
                  ),
                  tabPanel("Significant Difference Gene Venn Diagram",
                    plotOutput("venn_plot_2"),
                    downloadButton("download_venn_plot_2", "Download Venn Diagram 2")
                  ),
                  tabPanel("Significant Difference Gene Venn Diagram (Gradual Change)",
                    plotOutput("venn_plot_3"),
                    downloadButton("download_venn_plot_3", "Download Venn Diagram 3")
                  ),
                  tabPanel("Significant Difference Gene Heatmap",
                    plotOutput("significant_heatmap"),
                    downloadButton("download_significant_heatmap", "Download Heatmap of Significant Differences")
                  ),
                  tabPanel("CRLM Significant Difference Gene Heatmap",
                    plotOutput("crlm_significant_heatmap"),
                    downloadButton("download_crlm_significant_heatmap", "Download Heatmap of CRLM Significant Differences")
                  ),
                  tabPanel("CRLM Significant Difference Gene Heatmap (Sorted by Glycan Type)",
                    plotOutput("crlm_sorted_heatmap"),
                    downloadButton("download_crlm_sorted_heatmap", "Download Sorted Heatmap of CRLM Significant Differences")
                  ),
                  tabPanel("Gene List",
                    DTOutput("venn_genes_table"),
                    downloadButton("download_venn_genes", "Download Gene List")
                  ),
                  tabPanel("Combined Venn Diagram Display",
                    tags$p("Below are three Venn Diagram combinations, the heatmap of significant differences is shown separately below."),
                    plotOutput("venn_combined_plot"),
                    downloadButton("download_venn_combined", "Download Combined Venn Diagram")
                  )
                )
              )
            )
          ),
          tabPanel("LDA Plot and GO Analysis",
            fluidRow(
              box(width = 12, title = "LDA Plot and GO Analysis Input",
                actionButton("run_lda_go_analysis", "Run LDA Plot and GO Analysis"),
                helpText('This analysis will use the glycan-glycantype.csv file uploaded on the "Data Upload" page, as well as the two files in the "Data Upload" tab.')
              )
            ),
            fluidRow(
              box(width = 12, title = "LDA Plot and GO Analysis Results",
                tabsetPanel(
                  tabPanel("PCA Analysis",
                    plotOutput("pca_plot"),
                    downloadButton("download_pca_plot", "Download PCA Plot")
                  ),
                  tabPanel("LDA Analysis",
                    plotOutput("lda_plot"),
                    downloadButton("download_lda_plot", "Download LDA Plot")
                  ),
                  tabPanel("SVM Analysis",
                    plotOutput("svm_plot"),
                    downloadButton("download_svm_plot", "Download SVM Plot")
                  ),
                  tabPanel("GO Enrichment Analysis (Bar Plot)",
                    plotOutput("go_barplot"),
                    downloadButton("download_go_barplot", "Download GO Bar Plot")
                  ),
                  tabPanel("GO Enrichment Analysis (Dot Plot)",
                    plotOutput("go_dotplot"),
                    downloadButton("download_go_dotplot", "Download GO Dot Plot")
                  ),
                  tabPanel("Glycan Type Distribution (Pie Chart)",
                    plotOutput("glycan_pie_plot"),
                    downloadButton("download_glycan_pie_plot", "Download Glycan Type Pie Chart")
                  ),
                  tabPanel("Glycan Type Distribution (Bar Chart)",
                    plotOutput("glycan_bar_plot"),
                    downloadButton("download_glycan_bar_plot", "Download Glycan Type Bar Chart")
                  ),
                  tabPanel("GO Enrichment Results",
                    DTOutput("go_enrichment_table"),
                    downloadButton("download_go_enrichment", "Download GO Enrichment Results")
                  ),
                  tabPanel("Glycan Type Statistics",
                    DTOutput("glycan_type_table"),
                    downloadButton("download_glycan_type", "Download Glycan Type Statistics")
                  )
                )
              )
            )
          ),
          tabPanel("KEGG Analysis",
            fluidRow(
              box(width = 12, title = "KEGG Analysis Input",
                actionButton("run_kegg_analysis", "Run KEGG Analysis"),
                helpText('This analysis will use the glycan-glycantype.csv file uploaded on the "Data Upload" page, as well as the two files in the "Data Upload" tab.')
              )
            ),
            fluidRow(
              box(width = 12, title = "KEGG Analysis Results",
                tabsetPanel(
                  tabPanel("KEGG Table Combination Display",
                    plotOutput("kegg_combined_plot"),
                    downloadButton("download_kegg_combined", "Download Combination Table")
                  ),
                  tabPanel("KEGG Enrichment Barplot",
                    plotOutput("kegg_barplot"),
                    downloadButton("download_kegg_barplot", "Download Barplot")
                  ),
                  tabPanel("KEGG emapplot",
                    plotOutput("kegg_emapplot"),
                    downloadButton("download_kegg_emapplot", "Download emapplot")
                  ),
                  tabPanel("KEGG cnetplot",
                    plotOutput("kegg_cnetplot"),
                    downloadButton("download_kegg_cnetplot", "Download cnetplot")
                  ),
                  tabPanel("KEGG Enrichment Results Table",
                    DTOutput("kegg_table"),
                    downloadButton("download_kegg_table", "Download KEGG Enrichment Table")
                  )
                )
              )
            )
          ),
          tabPanel("Volcano Plot Analysis",
            fluidRow(
              box(width = 12, title = "Volcano Plot Analysis Input",
                actionButton("run_volcano_analysis", "Run Volcano Plot Analysis"),
                helpText('This analysis will use the glycan-glycantype.csv file uploaded on the "Data Upload" page, as well as the two files in the "Data Upload" tab.')
              )
            ),
            fluidRow(
              box(width = 12, title = "Volcano Plot Analysis Results",
                tabsetPanel(
                  tabPanel("Volcano Plot Combination Display",
                    plotOutput("volcano_combined_plot"),
                    downloadButton("download_volcano_combined", "Download Combined Volcano Plot")
                  ),
                  tabPanel("CRLM vs CRC Volcano Plot (All Marked)",
                    plotOutput("volcano_crc_1"),
                    downloadButton("download_volcano_crc_1", "Download CRLM vs CRC Volcano Plot 1")
                  ),
                  tabPanel("CRLM vs CRC Volcano Plot (High Expression Marked)",
                    plotOutput("volcano_crc_2"),
                    downloadButton("download_volcano_crc_2", "Download CRLM vs CRC Volcano Plot 2")
                  ),
                  tabPanel("CRLM vs LIHC Volcano Plot (All Marked)",
                    plotOutput("volcano_lihc_1"),
                    downloadButton("download_volcano_lihc_1", "Download CRLM vs LIHC Volcano Plot 1")
                  ),
                  tabPanel("CRLM vs LIHC Volcano Plot (High Expression Marked)",
                    plotOutput("volcano_lihc_2"),
                    downloadButton("download_volcano_lihc_2", "Download CRLM vs LIHC Volcano Plot 2")
                  ),
                  tabPanel("Protein Expression Analysis",
                    selectInput("protein_id", "Select Protein ID:", 
                               choices = c("P02763", "P01009", "P01008", "P05155", "P12109", "P02461", 
                                          "Q99715", "P05997", "Q05707", "P12110", "P12111", "P06731", 
                                          "P13688", "P07339")),
                    selectInput("plot_type", "Select Chart Type:", choices = c("violin" = "violin", "box" = "box")),
                    plotOutput("protein_expression_plot"),
                    downloadButton("download_protein_plot", "Download Protein Expression Plot")
                  ),
                  tabPanel("Volcano Plot Data",
                    DTOutput("volcano_data_table"),
                    downloadButton("download_volcano_data", "Download Volcano Plot Data")
                  ),
                  tabPanel("Cluster Data",
                    DTOutput("cluster_data_table"),
                    downloadButton("download_cluster_data", "Download Cluster Data")
                  )
                )
              )
            )
          )
        )
      ),
      
      # Results export page
      tabItem(tabName = "export",
              fluidRow(
                box(width = 12,
                    title = "Export Options",
                    checkboxGroupInput("export_items",
                                       "Select Content to Export:",
                                       choices = c(
                                         "Glycomics Analysis Results" = "glycomics",
                                         "Glycopeptide Analysis Results" = "glycopeptide"
                                       )),
                    downloadButton("download_all", "Export Selected Results")
                )
              )
      ),
      
      # Help page
      tabItem(tabName = "help",
              fluidRow(
                box(width = 12, title = "Help",
                    tabsetPanel(
                      tabPanel("Usage Flow",
                               h3("Step 1: Data Upload"),
                               p("In the \"Data Upload\" page, upload three required CSV files in sequence:"),
                               tags$ul(
                                 tags$li(strong("Group Normalization Data File: "), "contains glycopeptide identification and quantification information."),
                                 tags$li(strong("Patient Information File: "), "contains the correspondence between samples and patients."),
                                 tags$li(strong("Glycan Type File: "), "contains the correspondence between glycan composition and glycan type.")
                               ),
                               h3("Step 2: Data Preview and Quality Control"),
                               p("After uploading data, switch to the \"Data Preview and Quality Control\" page to check:"),
                               tags$ul(
                                 tags$li("Whether the data table is correctly displayed."),
                                 tags$li("Sample matching status, ensuring all samples are correctly matched."),
                                 tags$li("Data distribution graph, understanding overall data situation.")
                               ),
                               h3("Step 3: Execute Analysis"),
                               p("After confirming the data is correct, you can choose to perform \"Glycomics Analysis\", \"Glycoproteomics Qualitative Analysis\", or \"Glycoproteomics Quantitative Analysis\"."),
                               tags$ul(
                                 tags$li(strong("Glycomics Analysis: "), "This analysis module focuses on the change of glycan (Glycan) as a whole."),
                                 tags$li(strong("Glycoproteomics Qualitative Analysis: "), "This analysis module focuses on specific glycopeptide (Glycopeptide) molecules."),
                                 tags$li(
                                   strong("Glycoproteomics Quantitative Analysis: "),
                                   "This module includes 8 submodules, supporting data upload from expression matrix and group information, heatmap analysis, regional heatmap, Venn diagram, LDA/GO, KEGG, ROC, and multi-dimensional quantitative analysis."
                                 ),
                                 tags$ul(
                                   tags$li("Data Upload: Upload Imputed Missing Values _ zscore After.csv and sample_group_PCA.csv, supporting data preview."),
                                   tags$li("Heatmap Analysis: Various clustering and standardized heatmaps, showing expression patterns among different groups."),
                                   tags$li("Regional Heatmap Analysis: Regional expression differences and z-score heatmaps."),
                                   tags$li("Venn Diagram Analysis: Gene overlap and heatmaps among various screening methods."),
                                   tags$li("LDA/GO Analysis: PCA, LDA, SVM dimensionality reduction and GO enrichment."),
                                   tags$li("KEGG Analysis: Pathway enrichment, network diagrams, and tables."),
                                   tags$li("Volcano Plot Analysis: CRLM vs CRC/LIHC differences, protein expression, data tables.")
                                 )
                               ),
                               h3("Step 4: View and Export Results"),
                               p("After analysis is complete, view the generated charts and data on the corresponding page. Finally, you can package and download all the analysis results you need on the \"Export Results\" page.")
                      ),
                      tabPanel("Chart Explanation",
                               h4("Data Preview and Quality Control"),
                               tags$ul(
                                 tags$li(strong("Intensity Distribution: "), "shows the distribution of `final_scale_intensity` across all samples, which can be used to determine whether further conversion or standardization is needed."),
                                 tags$li(strong("Missing Value Overview: "), "shows the proportion of missing values in each sample, helping to identify low-quality samples.")
                               ),
                               h4("Glycomics Analysis"),
                               tags$ul(
                                 tags$li(strong("Molar Mass Distribution: "), "shows the distribution of glycan based on their molar mass across different groups. Y-axis represents the signal intensity after normalization."),
                                 tags$li(strong("Molar Mass Distribution (with Significant Annotation): "), "on the basis of the previous figure, glycan with significant differences in inter-group comparison (p < 0.05) are marked with red text."),
                                 tags$li(strong("Significant Difference Volcano Plot: "), "a volcano plot showing significant differences between groups, X-axis is the fold change (log2 Fold Change), Y-axis is -log10 of p-value. The points closer to the top are more significant.")
                               ),
                               h4("Glycoproteomics Qualitative Analysis"),
                               tags$ul(
                                 tags$li(strong("Sample Distribution Boxplot: "), "shows the distribution of sample intensity, which can be used to assess whether the overall intensity among samples is consistent."),
                                 tags$li(strong("Correlation Heatmap: "), "shows the correlation matrix between samples. The darker the color, the higher the similarity between samples. It can be used to assess the consistency of repeated samples."),
                                 tags$li(strong("PCA Analysis Plot: "), "the score plot of principal component analysis (PCA), used for sample clustering. Points closer together represent more similar samples."),
                                 tags$li(strong("Glycosylation Site Distribution: "), "counts the number of different glycan types connected on each peptide."),
                                 tags$li(strong("Sankey Diagram: "), "shows the flow of cancer type, glycan type, and number of glycosylation sites among them."),
                                 tags$li(strong("Unique Counts: "), "counts the number of unique proteins, glycans, peptides, and glycopeptides found in different cancer types.")
                               )
                      ),
                      tabPanel("Data Format Requirements", includeMarkdown("data_format_examples.md"))
                    )
                )
              )
      )
    )
  )
)

# Define server logic
server <- function(input, output, session) {
  # Reactive values
  zujian_data <- reactiveVal(NULL)
  patient_info_data <- reactiveVal(NULL)
  glycan_type_data <- reactiveVal(NULL)
  glycomics_results <- reactiveVal(NULL)
  glycopeptide_results <- reactiveVal(NULL)
  igp_expr_data <- reactiveVal(NULL)
  igp_group_data <- reactiveVal(NULL)
  igp_heatmap_result <- reactiveVal(NULL)
  regional_analysis_results <- reactiveVal(NULL)
  venn_analysis_results <- reactiveVal(NULL)
  lda_go_analysis_results <- reactiveVal(NULL)
  kegg_analysis_results <- reactiveVal(NULL)
  volcano_analysis_results <- reactiveVal(NULL)
  
  # Home page navigation button events
  observeEvent(input$go_data_processing, {
    updateTabItems(session, "tabs", "data_processing")
  })
  
  observeEvent(input$go_data_analysis, {
    updateTabItems(session, "tabs", "data_upload")
  })
  
  # Load data after upload
  observeEvent(input$zujian_data, {
    req(input$zujian_data)
    data <- read.csv(input$zujian_data$datapath, header = TRUE)
    zujian_data(data)
    showNotification("Group Normalization Data Loaded Successfully", type = "message")
  })
  
  observeEvent(input$patient_info, {
    req(input$patient_info)
    data <- read.csv(input$patient_info$datapath, header = TRUE)
    patient_info_data(data)
    showNotification("Patient Information Data Loaded Successfully", type = "message")
  })
  
  observeEvent(input$glycan_type, {
    req(input$glycan_type)
    data <- read.csv(input$glycan_type$datapath, header = TRUE)
    glycan_type_data(data)
    showNotification("Glycan Type Data Loaded Successfully", type = "message")
  })
  
  # Data upload status
  output$upload_status <- renderText({
    status <- c()
    if (!is.null(zujian_data())) {
      status <- c(status, paste("✓ Group Normalization Data Uploaded: ", nrow(zujian_data()), " rows x", ncol(zujian_data()), " columns"))
    } else {
      status <- c(status, "✗ Group Normalization Data: Not Uploaded")
    }
    
    if (!is.null(patient_info_data())) {
      status <- c(status, paste("✓ Patient Information Data Uploaded: ", nrow(patient_info_data()), " rows x", ncol(patient_info_data()), " columns"))
    } else {
      status <- c(status, "✗ Patient Information Data: Not Uploaded")
    }
    
    if (!is.null(glycan_type_data())) {
      status <- c(status, paste("✓ Glycan Type Data Uploaded: ", nrow(glycan_type_data()), " rows x", ncol(glycan_type_data()), " columns"))
    } else {
      status <- c(status, "✗ Glycan Type Data: Not Uploaded")
    }
    
    paste(status, collapse = "\n")
  })
  
  # Data preview
  output$zujian_data_preview <- renderDT({
    req(zujian_data())
    datatable(head(zujian_data(), 100), options = list(scrollX = TRUE, pageLength = 5))
  })
  
  output$patient_info_preview <- renderDT({
    req(patient_info_data())
    datatable(patient_info_data(), options = list(scrollX = TRUE, pageLength = 5))
  })
  
  output$glycan_type_preview <- renderDT({
    req(glycan_type_data())
    datatable(glycan_type_data(), options = list(scrollX = TRUE, pageLength = 5))
  })
  
  # Data quality control and visualization
  output$sample_matching_status <- renderText({
    req(zujian_data(), patient_info_data())
    norm_ids <- unique(zujian_data()$sample_id)
    patient_ids <- unique(patient_info_data()$Sample)
    
    report <- c("=== Sample Matching Status ===")
    report <- c(report, paste("Number of Independent Samples in Group Normalization Data:", length(norm_ids)))
    report <- c(report, paste("Number of Independent Samples in Patient Information Data:", length(patient_ids)))
    
    unmatched_in_patient <- setdiff(patient_ids, norm_ids)
    if (length(unmatched_in_patient) > 0) {
      report <- c(report, paste("✗ Warning: ", length(unmatched_in_patient), " samples in patient information cannot be found in group normalization data."))
    }
    
    unmatched_in_norm <- setdiff(norm_ids, patient_ids)
    if (length(unmatched_in_norm) > 0) {
      report <- c(report, paste("✗ Warning: ", length(unmatched_in_norm), " samples in group normalization data cannot be found in patient information."))
    }
    
    matched_count <- length(intersect(norm_ids, patient_ids))
    report <- c(report, paste("✓ Successful Matches:", matched_count))
    
    if (matched_count > 0) {
      report <- c(report, "✓ Sample Matching Successful, Analysis Can Proceed.")
    } else {
      report <- c(report, "✗ Error: No Samples Successfully Matched, Please Check Whether Sample IDs Are Consistent.")
    }
    paste(report, collapse = "\n")
  })
  
  output$intensity_distribution_plot <- renderPlot({
    req(zujian_data())
    ggplot(zujian_data(), aes(x = final_scale_intensity)) +
      geom_histogram(bins = 50, fill = "lightblue", color = "black") +
      labs(title = "Data Intensity Distribution (final_scale_intensity)", 
           x = "Log2(Intensity)", 
           y = "Frequency") +
      theme_classic()
  })
  
  output$missing_value_plot <- renderPlot({
    req(zujian_data())
    missing_data <- zujian_data() %>%
      group_by(sample_id) %>%
      summarise(missing_ratio = sum(is.na(final_scale_intensity) | final_scale_intensity == 0) / n(), .groups = 'drop')
    
    ggplot(missing_data, aes(x = reorder(sample_id, -missing_ratio), y = missing_ratio)) +
      geom_bar(stat = "identity", fill = "steelblue") +
      scale_y_continuous(labels = scales::percent) +
      labs(title = "Proportion of Missing Values per Sample",
           x = "Sample ID",
           y = "Missing Value Proportion") +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 8))
  })
  
  # Glycomics analysis
  observeEvent(input$run_glycomics, {
    req(zujian_data(), patient_info_data())
    withProgress(message = 'Running Glycomics Analysis...', {
      # Check sample ID matching
      norm_ids <- unique(zujian_data()$sample_id)
      patient_ids <- unique(patient_info_data()$Sample)
      not_in_norm <- setdiff(patient_ids, norm_ids)
      not_in_patient <- setdiff(norm_ids, patient_ids)
      if (length(not_in_norm) > 0) {
        showNotification(paste0("Samples below are not found in group normalization data: ", paste(not_in_norm, collapse=", ")), type = "warning")
      }
      if (length(not_in_patient) > 0) {
        showNotification(paste0("Samples below are not found in patient information table: ", paste(not_in_patient, collapse=", ")), type = "warning")
      }
      
      # Check CRC grouping logic
      crc_in_patient <- any(grepl("CRC", patient_info_data()$patients_repeat))
      norm_patient_info <- merge(
        data.frame(sample_id = norm_ids),
        patient_info_data(),
        by.x = "sample_id", by.y = "Sample",
        all.x = TRUE
      )
      crc_in_norm <- any(grepl("CRC", norm_patient_info$patients_repeat))
      if (!crc_in_norm || !crc_in_patient) {
        showNotification("Warning: CRC grouping is missing in group normalization data or patient information table, CRC will not be displayed in the figure!", type = "error")
      }
      
      glyco_res <- glycomics_analysis_20250521(zujian_data(), patient_info_data())
      
      # Ensure all cancertype groups participate in facet_wrap
      plot_data <- glyco_res$normalized
      plot_data$cancertype <- factor(plot_data$cancertype, levels = glyco_res$cancertype_levels)
      
      output$glycomics_mass_plot <- renderPlot({
        x_ticks <- sort(unique(plot_data$Glyco_Molar_Mass))
        ggplot(plot_data, aes(x = Glyco_Molar_Mass, y = Mean_Value_normalized)) +
          geom_segment(aes(xend = Glyco_Molar_Mass, yend = 0), size = 1) +
          facet_wrap(~ cancertype, scales = "fixed", ncol = 1, drop = FALSE) +
          scale_x_continuous(limits = c(570, 6500), breaks = x_ticks) +
          labs(x = "Glyco_Molar_Mass", y = "Mean_Value") +
          theme_minimal() +
          theme(
            strip.text = element_text(size = 12),
            axis.text = element_text(size = 10),
            axis.title = element_text(size = 12)
          )
      })
      
      # Render molar mass distribution (with significance annotation)
      output$glycomics_mass_sig_plot <- renderPlot({
        sig <- glyco_res$significant
        x_ticks <- sort(unique(plot_data$Glyco_Molar_Mass))
        p <- ggplot(plot_data, aes(x = Glyco_Molar_Mass, y = Mean_Value_normalized)) +
          geom_segment(aes(xend = Glyco_Molar_Mass, yend = 0), size = 1) +
          facet_wrap(~ cancertype, scales = "fixed", ncol = 1, drop = FALSE) +
          scale_x_continuous(limits = c(570, 5000), breaks = x_ticks) +
          labs(x = "Glyco_Molar_Mass", y = "Mean_Value") +
          theme_minimal() +
          theme(
            strip.text = element_text(size = 12),
            axis.text = element_text(size = 10),
            axis.title = element_text(size = 12)
          )
        if (nrow(sig) > 0) {
          p <- p + geom_text(data = subset(plot_data, Glyco_Molar_Mass %in% sig$Glyco_Molar_Mass),
                             aes(label = Glyco_Molar_Mass), color = "red", size = 3, vjust = -0.5)
        }
        p
      })
      
      # Render significant differences
      output$glycomics_significant_plot <- renderPlot({
        sig <- glyco_res$significant
        if (nrow(sig) > 0) {
          ggplot(sig, aes(x = Glyco_Molar_Mass, y = -log10(p_value), color = comparison)) +
            geom_point() +
            labs(x = "Glyco_Molar_Mass", y = "-log10(p_value)", color = "Comparison") +
            theme_minimal()
        } else {
          plot(1, 1, type = "n", axes = FALSE, xlab = "", ylab = "")
          text(1, 1, "No Significant Differences Found", cex = 1.5)
        }
      })
      
      # Render t-test results table
      output$glycomics_ttest_table <- renderDT({
        datatable(glyco_res$significant, options = list(scrollX = TRUE))
      })
      
      # Render data table
      output$glycomics_table <- renderDT({
        datatable(glyco_res$normalized, options = list(scrollX = TRUE))
      })
      
      glycomics_results(glyco_res)
    })
  })
  
  # Glycopeptide analysis
  observeEvent(input$run_glycopeptide, {
    req(zujian_data(), patient_info_data(), glycan_type_data())
    withProgress(message = 'Running Glycopeptide Analysis...', {
      glycopep_res <- glycopeptide_analysis_20250521(zujian_data(), patient_info_data(), glycan_type_data())
      
      # Unique statistics bar chart
      output$unique_counts_plot <- renderPlot({
        ggplot(glycopep_res$unique_counts, aes(x = factor(cancer_type, levels = glycopep_res$cancer_type_levels), y = count, fill = unique_type)) +
          geom_bar(stat = "identity", position = "dodge") +
          labs(title = "Unique Counts by Cancer Type",
               x = "Cancer Type",
               y = "Count",
               fill = "Unique Type") +
          geom_text(aes(label = count),
                    position = position_dodge(width = 0.9),
                    vjust = -0.3) +
          theme_classic()
      })
      
      # Sample distribution boxplot
      output$glycopeptide_boxplot <- renderPlot({
        ggplot(glycopep_res$long, aes(x = Sample, y = Value)) +
          geom_boxplot() +
          theme_classic() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
          labs(title = "Glycopeptide Data Distribution", x = "Sample", y = "Intensity Value")
      })
      
      # Peptide glycan type histogram
      output$glycan_hist <- renderPlot({
        ggplot(glycopep_res$glycan_counts_per_peptide, aes(x = unique_glycans)) +
          geom_histogram(binwidth = 1, fill = "blue", color = "black", alpha = 0.7) +
          labs(title = "Distribution of Unique Glycan Types per Peptide",
               x = "Number of Unique Glycan Types",
               y = "Frequency") +
          theme_classic()
      })
      
      # Peptide glycan type density plot
      output$glycan_density <- renderPlot({
        ggplot(glycopep_res$glycan_counts_per_peptide, aes(x = unique_glycans)) +
          geom_density(fill = "blue", alpha = 0.7) +
          labs(title = "Density of Unique Glycan Types per Peptide",
               x = "Number of Unique Glycan Types",
               y = "Density") +
          theme_classic()
      })
      
      # Glycopeptide type bar chart
      output$glycan_type_bar <- renderPlot({
        ggplot(glycopep_res$glycan_type_counts, aes(x = glycan.type, y = glycopeptides_count)) +
          geom_bar(stat = "identity", fill = "blue", color = "black", alpha = 0.7) +
          labs(title = "Number of Glycopeptides per Glycan Type",
               x = "Glycan Type",
               y = "Number of Glycopeptides") +
          theme_classic() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1))
      })
      
      # Glycosylation site distribution
      output$glycosite_distribution <- renderPlot({
        ggplot(glycopep_res$glycan_counts_per_peptide, aes(x = unique_glycans)) +
          geom_histogram(binwidth = 1, fill = "green", color = "black", alpha = 0.7) +
          labs(title = "Distribution of Glycosylation Sites per Peptide",
               x = "Number of Glycosylation Sites",
               y = "Frequency") +
          theme_classic()
      })
      
      # Sankey diagram
      output$sankey_plot <- renderPlot({
        library(ggalluvial)
        df <- glycopep_res$df_for_dingxing_analysis_long
        glyco_glycan_sites <- glycopep_res$glyco_glycan_sites
        # Merge site information
        df_sankey <- merge(df, glyco_glycan_sites, by = c("peptides", "glycan", "glycan.type"))
        # Generate site labels
        df_sankey <- df_sankey %>%
          dplyr::mutate(sites_label = dplyr::case_when(
            sites == 1 ~ "peptides contain 1 glycan composition",
            sites == 2 ~ "peptides contain 2 glycan composition",
            sites == 3 ~ "peptides contain 3 glycan composition",
            sites == 4 ~ "peptides contain 4 glycan composition",
            sites == 5 ~ "peptides contain 5 glycan composition",
            sites == 6 ~ "peptides contain 6 glycan composition",
            sites > 6 ~ "peptides contain >6 glycan composition"
          ))
        # Statistical counting
        df_sankey_count <- df_sankey %>%
          dplyr::group_by(cancer_type, glycan.type, sites_label) %>%
          dplyr::summarise(count = n(), .groups = "drop")
        # Plotting
        ggplot(df_sankey_count,
               aes(axis1 = cancer_type, axis2 = glycan.type, axis3 = sites_label, y = count)) +
          geom_alluvium(aes(fill = cancer_type), width = 1/12) +
          geom_stratum(width = 1/8, fill = "white", color = "black") +
          geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 4) +
          scale_x_discrete(limits = c("Cancer Type", "Glycan Type", "Sites")) +
          scale_fill_manual(values = c("CRC" = "#F8766D", "CRLM" = "#00BA38", "LIHC" = "#619CFF")) +
          theme_minimal() +
          labs(title = "Sankey Diagram of Glycan Types, Cancer Types, and Sites",
               x = "Categories", y = "Count", fill = "cancer_type")
      })
      
      # PCA analysis
      output$pca_plot <- renderPlot({
        library(ggfortify)
        df_wide <- glycopep_res$wide
        # Keep only numeric sample columns
        df_wide <- df_wide %>% dplyr::select(Glycopeptides, where(is.numeric))
        # Filter rows with high missing ratio
        df_wide <- df_wide[rowSums(is.na(df_wide[,-1]))/ncol(df_wide[,-1]) <= 0.5, ]
        # KNN imputation
        df_filled <- VIM::kNN(df_wide, variable = setdiff(colnames(df_wide), "Glycopeptides"))
        df_filled <- dplyr::select(df_filled, -ends_with("_imp"))
        data_for_pca <- df_filled
        rownames(data_for_pca) <- data_for_pca$Glycopeptides
        data_for_pca <- data_for_pca %>% dplyr::select(-Glycopeptides)
        data_for_pca <- data_for_pca[, order(names(data_for_pca))]
        data_for_pca <- data_for_pca[, apply(data_for_pca, 2, function(x) var(x, na.rm = TRUE) > 0)]
        data_for_pca <- data_for_pca[complete.cases(data_for_pca), ]
        df_scaled <- scale(data_for_pca)
        pca_result <- prcomp(df_scaled, center = TRUE, scale. = TRUE)
        autoplot(pca_result, main = "PCA Analysis") +
          theme_classic() +
          labs(x = "Principal Component 1 Loadings (PC1)", y = "Principal Component 2 Loadings (PC2)")
      })
      
      # Correlation heatmap
      output$correlation_heatmap <- renderPlot({
        library(ggcorrplot)
        df_wide <- glycopep_res$wide
        df_wide[df_wide == "NA"] <- NA
        df_wide <- df_wide %>% dplyr::select(Glycopeptides, where(is.numeric))
        df_wide <- df_wide[rowSums(is.na(df_wide[,-1]))/ncol(df_wide[,-1]) <= 0.5, ]
        df_filled <- VIM::kNN(df_wide, variable = setdiff(colnames(df_wide), "Glycopeptides"))
        df_filled <- dplyr::select(df_filled, -ends_with("_imp"))
        data_for_cor <- df_filled
        rownames(data_for_cor) <- data_for_cor$Glycopeptides
        data_for_cor <- data_for_cor %>% dplyr::select(-Glycopeptides)
        data_for_cor <- data_for_cor[, order(names(data_for_cor))]
        data_for_cor <- data_for_cor[, apply(data_for_cor, 2, function(x) var(x, na.rm = TRUE) > 0)]
        data_for_cor <- data_for_cor[complete.cases(data_for_cor), ]
        correlation_matrix <- cor(data_for_cor, use = "pairwise.complete.obs")
        ggcorrplot(correlation_matrix, 
                   method = "square", 
                   type = "lower", 
                   lab = FALSE, 
                   title = "Correlation Matrix Heatmap", 
                   ggtheme = theme_classic())
      })
      
      # Data table
      output$glycopeptide_table <- renderDT({
        datatable(glycopep_res$wide, options = list(scrollX = TRUE))
      })
      
      # New: Three intermediate result tables
      output$knn_input_table <- renderDT({
        datatable(glycopep_res$knn_input, options = list(scrollX = TRUE))
      })
      
      output$correlation_matrix_table <- renderDT({
        datatable(as.data.frame(glycopep_res$correlation_matrix), options = list(scrollX = TRUE))
      })
      
      output$zscore_table <- renderDT({
        datatable(glycopep_res$zscore, options = list(scrollX = TRUE))
      })
      
      glycopeptide_results(glycopep_res)
    })
  })
  
  # Glycopeptide quantitative analysis data upload
  observeEvent(input$igp_expr_file, {
    req(input$igp_expr_file)
    data <- read.csv(input$igp_expr_file$datapath, header = FALSE, stringsAsFactors = FALSE)
    igp_expr_data(data)
    showNotification("Expression Data Loaded Successfully", type = "message")
  })
  
  observeEvent(input$igp_group_file, {
    req(input$igp_group_file)
    data <- read.csv(input$igp_group_file$datapath, header = TRUE, stringsAsFactors = FALSE)
    igp_group_data(data)
    showNotification("Group Information Loaded Successfully", type = "message")
  })
  
  output$igp_expr_preview <- renderDT({ req(igp_expr_data()); datatable(head(igp_expr_data(), 10), options = list(scrollX = TRUE)) })
  output$igp_group_preview <- renderDT({ req(igp_group_data()); datatable(head(igp_group_data(), 10), options = list(scrollX = TRUE)) })
  
  # Glycopeptide quantitative analysis heatmap
  observeEvent(input$run_igp_heatmap, {
    req(igp_expr_data(), igp_group_data())
    withProgress(message = 'Analyzing and Plotting Heatmap...', {
      res <- igp_quant_heatmap_analysis(igp_expr_data(), igp_group_data())
      igp_heatmap_result(res)
    })
  })
  
  output$igp_heatmap_crc <- renderPlot({
    req(igp_heatmap_result())
    if (!is.null(igp_heatmap_result()$crc_heatmap)) {
      # If group annotation is needed, use group_levels
      igp_heatmap_result()$crc_heatmap
    }
  })
  
  output$igp_heatmap_lihc <- renderPlot({
    req(igp_heatmap_result())
    if (!is.null(igp_heatmap_result()$lihc_heatmap)) {
      igp_heatmap_result()$lihc_heatmap
    }
  })
  
  output$igp_heatmap_combined <- renderPlot({
    req(igp_heatmap_result())
    igp_heatmap_result()$combined_heatmap
  })
  
  output$igp_heatmap_zscore <- renderPlot({
    req(igp_heatmap_result())
    igp_heatmap_result()$zscore_heatmap
  })
  
  output$igp_heatmap_combined_nocluster <- renderPlot({
    req(igp_heatmap_result())
    igp_heatmap_result()$combined_heatmap_nocluster
  })
  
  output$igp_heatmap_zscore_cluster <- renderPlot({
    req(igp_heatmap_result())
    igp_heatmap_result()$zscore_heatmap_cluster
  })
  
  output$igp_heatmap_zscore_nocluster <- renderPlot({
    req(igp_heatmap_result())
    igp_heatmap_result()$zscore_heatmap_nocluster
  })
  
  output$download_igp_heatmap_crc <- downloadHandler(
    filename = function() { "IGP_heatmap_crc.png" },
    content = function(file) {
      png(file, width = 1200, height = 900, res = 120)
      res <- igp_heatmap_result()
      if (!is.null(res) && !is.null(res$crc_heatmap)) print(res$crc_heatmap)
      dev.off()
    }
  )
  
  output$download_igp_heatmap_lihc <- downloadHandler(
    filename = function() { "IGP_heatmap_lihc.png" },
    content = function(file) {
      png(file, width = 1200, height = 900, res = 120)
      res <- igp_heatmap_result()
      if (!is.null(res) && !is.null(res$lihc_heatmap)) print(res$lihc_heatmap)
      dev.off()
    }
  )
  
  output$download_igp_heatmap_combined <- downloadHandler(
    filename = function() { "IGP_heatmap_combined.png" },
    content = function(file) {
      png(file, width = 1200, height = 900, res = 120)
      res <- igp_heatmap_result()
      if (!is.null(res)) print(res$combined_heatmap)
      dev.off()
    }
  )
  
  output$download_igp_heatmap_zscore <- downloadHandler(
    filename = function() { "IGP_heatmap_zscore.png" },
    content = function(file) {
      png(file, width = 1200, height = 900, res = 120)
      res <- igp_heatmap_result()
      if (!is.null(res)) print(res$zscore_heatmap)
      dev.off()
    }
  )
  
  output$download_igp_heatmap_combined_nocluster <- downloadHandler(
    filename = function() { "IGP_heatmap_combined_nocluster.png" },
    content = function(file) {
      png(file, width = 1200, height = 900, res = 120)
      res <- igp_heatmap_result()
      if (!is.null(res)) print(res$combined_heatmap_nocluster)
      dev.off()
    }
  )
  
  output$download_igp_heatmap_zscore_cluster <- downloadHandler(
    filename = function() { "IGP_heatmap_zscore_cluster.png" },
    content = function(file) {
      png(file, width = 1200, height = 900, res = 120)
      res <- igp_heatmap_result()
      if (!is.null(res)) print(res$zscore_heatmap_cluster)
      dev.off()
    }
  )
  
  output$download_igp_heatmap_zscore_nocluster <- downloadHandler(
    filename = function() { "IGP_heatmap_zscore_nocluster.png" },
    content = function(file) {
      png(file, width = 1200, height = 900, res = 120)
      res <- igp_heatmap_result()
      if (!is.null(res)) print(res$zscore_heatmap_nocluster)
      dev.off()
    }
  )
  
  # Download handling
  # Glycomics analysis related download buttons
  output$download_glycomics_data <- downloadHandler(
    filename = function() { "glycomics_normalized.csv" },
    content = function(file) {
      write.csv(glycomics_results()$normalized, file, row.names = FALSE)
    }
  )
  
  output$download_mass_plot <- downloadHandler(
    filename = function() { "glycomics_mass_plot.png" },
    content = function(file) {
      png(file, width = 1200, height = 800, res = 120)
      x_ticks <- sort(unique(glycomics_results()$normalized$Glyco_Molar_Mass))
      print(
        ggplot(glycomics_results()$normalized, aes(x = Glyco_Molar_Mass, y = Mean_Value_normalized)) +
          geom_segment(aes(xend = Glyco_Molar_Mass, yend = 0), size = 1) +
          facet_wrap(~ cancertype, scales = "fixed", ncol = 1) +
          scale_x_continuous(limits = c(570, 6500), breaks = x_ticks) +
          labs(x = "Glyco_Molar_Mass", y = "Mean_Value") +
          theme_minimal() +
          theme(
            strip.text = element_text(size = 12),
            axis.text = element_text(size = 10),
            axis.title = element_text(size = 12)
          )
      )
      dev.off()
    }
  )
  
  output$download_significant_plot <- downloadHandler(
    filename = function() { "glycomics_significant_plot.png" },
    content = function(file) {
      png(file, width = 1200, height = 800, res = 120)
      sig <- glycomics_results()$significant
      if (nrow(sig) > 0) {
        print(
          ggplot(sig, aes(x = Glyco_Molar_Mass, y = -log10(p_value), color = comparison)) +
            geom_point() +
            labs(x = "Glyco_Molar_Mass", y = "-log10(p_value)", color = "Comparison") +
            theme_minimal()
        )
      } else {
        plot(1, 1, type = "n", axes = FALSE, xlab = "", ylab = "")
        text(1, 1, "No Significant Differences Found", cex = 1.5)
      }
      dev.off()
    }
  )
  
  # Glycopeptide analysis related download buttons
  output$download_boxplot <- downloadHandler(
    filename = function() { "glycopeptide_boxplot.png" },
    content = function(file) {
      png(file, width = 1200, height = 800, res = 120)
      print(
        ggplot(glycopeptide_results()$unique_counts, aes(x = cancer_type, y = count, fill = unique_type)) +
          geom_bar(stat = "identity", position = "dodge") +
          labs(title = "Unique Counts by Cancer Type",
               x = "Cancer Type",
               y = "Count",
               fill = "Unique Type") +
          geom_text(aes(label = count),
                    position = position_dodge(width = 0.9),
                    vjust = -0.3) +
          theme_classic()
      )
      dev.off()
    }
  )
  
  output$download_heatmap <- downloadHandler(
    filename = function() { "correlation_heatmap.png" },
    content = function(file) {
      png(file, width = 1200, height = 800, res = 120)
      library(ggcorrplot)
      df_wide <- glycopeptide_results()$wide
      df_wide[df_wide == "NA"] <- NA
      df_wide <- df_wide %>% dplyr::select(Glycopeptides, where(is.numeric))
      df_wide <- df_wide[rowSums(is.na(df_wide[,-1]))/ncol(df_wide[,-1]) <= 0.5, ]
      df_filled <- VIM::kNN(df_wide, variable = setdiff(colnames(df_wide), "Glycopeptides"))
      df_filled <- dplyr::select(df_filled, -ends_with("_imp"))
      data_for_cor <- df_filled
      rownames(data_for_cor) <- data_for_cor$Glycopeptides
      data_for_cor <- data_for_cor %>% dplyr::select(-Glycopeptides)
      data_for_cor <- data_for_cor[, order(names(data_for_cor))]
      data_for_cor <- data_for_cor[, apply(data_for_cor, 2, function(x) var(x, na.rm = TRUE) > 0)]
      data_for_cor <- data_for_cor[complete.cases(data_for_cor), ]
      correlation_matrix <- cor(data_for_cor, use = "pairwise.complete.obs")
      print(
        ggcorrplot(correlation_matrix, 
                   method = "square", 
                   type = "lower", 
                   lab = FALSE, 
                   title = "Correlation Matrix Heatmap", 
                   ggtheme = theme_classic())
      )
      dev.off()
    }
  )
  
  output$download_pca <- downloadHandler(
    filename = function() { "pca_plot.png" },
    content = function(file) {
      png(file, width = 1200, height = 800, res = 120)
      library(ggfortify)
      df_wide <- glycopeptide_results()$wide
      df_wide <- df_wide %>% dplyr::select(Glycopeptides, where(is.numeric))
      df_wide <- df_wide[rowSums(is.na(df_wide[,-1]))/ncol(df_wide[,-1]) <= 0.5, ]
      df_filled <- VIM::kNN(df_wide, variable = setdiff(colnames(df_wide), "Glycopeptides"))
      df_filled <- dplyr::select(df_filled, -ends_with("_imp"))
      data_for_pca <- df_filled
      rownames(data_for_pca) <- data_for_pca$Glycopeptides
      data_for_pca <- data_for_pca %>% dplyr::select(-Glycopeptides)
      data_for_pca <- data_for_pca[, order(names(data_for_pca))]
      data_for_pca <- data_for_pca[, apply(data_for_pca, 2, function(x) var(x, na.rm = TRUE) > 0)]
      data_for_pca <- data_for_pca[complete.cases(data_for_pca), ]
      df_scaled <- scale(data_for_pca)
      pca_result <- prcomp(df_scaled, center = TRUE, scale. = TRUE)
      print(autoplot(pca_result, main = "PCA Analysis") + theme_classic())
      dev.off()
    }
  )
  
  output$download_glycosite <- downloadHandler(
    filename = function() { "glycosite_distribution.png" },
    content = function(file) {
      png(file, width = 1200, height = 800, res = 120)
      print(
        ggplot(glycopeptide_results()$glycan_counts_per_peptide, aes(x = unique_glycans)) +
          geom_histogram(binwidth = 1, fill = "blue", color = "black", alpha = 0.7) +
          labs(title = "Distribution of Unique Glycan Types per Peptide",
               x = "Number of Unique Glycan Types",
               y = "Frequency") +
          theme_classic()
      )
      dev.off()
    }
  )
  
  output$download_sankey <- downloadHandler(
    filename = function() { "sankey_plot.png" },
    content = function(file) {
      png(file, width = 1200, height = 800, res = 120)
      library(ggalluvial)
      df <- glycopeptide_results()$df_for_dingxing_analysis_long
      glyco_glycan_sites <- glycopeptide_results()$glyco_glycan_sites
      df_sankey <- merge(df, glyco_glycan_sites)
      df_sankey <- df_sankey %>%
        dplyr::mutate(sites_label = dplyr::case_when(
          sites == 1 ~ "peptides contain 1 glycan composition",
          sites == 2 ~ "peptides contain 2 glycan composition",
          sites == 3 ~ "peptides contain 3 glycan composition",
          sites == 4 ~ "peptides contain 4 glycan composition",
          sites == 5 ~ "peptides contain 5 glycan composition",
          sites == 6 ~ "peptides contain 6 glycan composition",
          sites > 6 ~ "peptides contain >6 glycan composition"
        )) %>%
        dplyr::mutate(cancer_glycan = paste(cancer_type, glycan.type, sep = "-"),
                      cancer_glycan_sites = paste(cancer_glycan, sites_label, sep = "_"))
      df_sankey_count <- df_sankey %>%
        dplyr::group_by(cancer_glycan_sites) %>%
        dplyr::summarise(counts = n())
      df_sankey <- merge(df_sankey, df_sankey_count)
      df_sankey <- unique(df_sankey)
      print(
        ggplot(df_sankey,
               aes(axis1 = cancer_type, axis2 = glycan.type, axis3 = sites_label, y = counts)) +
          geom_alluvium(aes(fill = cancer_type)) +
          geom_stratum() +
          geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
          scale_x_discrete(limits = c("Cancer Type", "Glycan Type", "Sites")) +
          theme_minimal() +
          labs(title = "Sankey Diagram of Glycan Types, Cancer Types, and Sites",
               x = "Categories", y = "Counts")
      )
      dev.off()
    }
  )
  
  output$download_glycopeptide_data <- downloadHandler(
    filename = function() { "glycopeptide_analysis_data.csv" },
    content = function(file) {
      write.csv(glycopeptide_results()$unique_counts, file, row.names = FALSE)
    }
  )
  
  # t-test results download button
  output$download_glycomics_ttest <- downloadHandler(
    filename = function() { "glycomics_ttest_results.csv" },
    content = function(file) {
      write.csv(glycomics_results()$significant, file, row.names = FALSE)
    }
  )
  
  # Molar mass distribution (with significance annotation) download button
  output$download_mass_sig_plot <- downloadHandler(
    filename = function() { "glycomics_mass_sig_plot.png" },
    content = function(file) {
      plot_data <- glycomics_results()$normalized
      sig <- glycomics_results()$significant
      all_types <- c("CRC", "CRLM", "LIHC")
      plot_data$cancertype <- factor(plot_data$cancertype, levels = all_types)
      x_ticks <- sort(unique(plot_data$Glyco_Molar_Mass))
      png(file, width = 1200, height = 800, res = 120)
      p <- ggplot(plot_data, aes(x = Glyco_Molar_Mass, y = Mean_Value_normalized)) +
        geom_segment(aes(xend = Glyco_Molar_Mass, yend = 0), size = 1) +
        facet_wrap(~ cancertype, scales = "fixed", ncol = 1, drop = FALSE) +
        scale_x_continuous(limits = c(570, 5000), breaks = x_ticks) +
        labs(x = "Glyco_Molar_Mass", y = "Mean_Value") +
        theme_minimal() +
        theme(
          strip.text = element_text(size = 12),
          axis.text = element_text(size = 10),
          axis.title = element_text(size = 12)
        )
      if (nrow(sig) > 0) {
        p <- p + geom_text(data = subset(plot_data, Glyco_Molar_Mass %in% sig$Glyco_Molar_Mass),
                           aes(label = Glyco_Molar_Mass), color = "red", size = 3, vjust = -0.5)
      }
      print(p)
      dev.off()
    }
  )
  
  # Data export
  output$download_all <- downloadHandler(
    filename = function() {
      paste("glycopeptide_analysis_", Sys.Date(), ".zip", sep = "")
    },
    content = function(file) {
      # Create temporary directory
      temp_dir <- tempdir()
      
      # Export data based on selection
      if ("glycomics" %in% input$export_items) {
        write.csv(glycomics_results()$normalized, file.path(temp_dir, "glycomics_results.csv"), row.names = FALSE)
        write.csv(glycomics_results()$significant, file.path(temp_dir, "glycomics_ttest_results.csv"), row.names = FALSE)
      }
      if ("glycopeptide" %in% input$export_items) {
        write.csv(glycopeptide_results()$unique_counts, file.path(temp_dir, "glycopeptide_results.csv"), row.names = FALSE)
        write.csv(glycopeptide_results()$wide, file.path(temp_dir, "glycopeptide_wide_data.csv"), row.names = FALSE)
      }
      
      # Create ZIP file
      files <- list.files(temp_dir, full.names = TRUE)
      zip(file, files)
    }
  )
  
  # New: Three intermediate result download buttons
  output$download_knn_input <- downloadHandler(
    filename = function() { "50% Effective Values Retained _ knn Imputation Before.csv" },
    content = function(file) {
      write.csv(glycopeptide_results()$knn_input, file, row.names = FALSE)
    }
  )
  
  output$download_correlation_matrix <- downloadHandler(
    filename = function() { "Sample Repeats Correlation.csv" },
    content = function(file) {
      write.csv(glycopeptide_results()$correlation_matrix, file)
    }
  )
  
  output$download_zscore <- downloadHandler(
    filename = function() { "Imputation of Missing Values _ zscore After.csv" },
    content = function(file) {
      write.csv(glycopeptide_results()$zscore, file, row.names = FALSE)
    }
  )

  # Regional heatmap analysis
  observeEvent(input$run_regional_analysis, {
    req(igp_expr_data(), igp_group_data(), glycan_type_data())
    withProgress(message = 'Running Regional Analysis...', {
      res <- regional_analysis_and_plots(igp_expr_data(), igp_group_data(), glycan_type_data())
      regional_analysis_results(res)
    })
  })

  output$regional_scatter_plot <- renderPlot({
    req(regional_analysis_results())
    regional_analysis_results()$scatter_plot
  })

  output$regional_counts_table <- renderDT({
    req(regional_analysis_results())
    datatable(regional_analysis_results()$region_counts, options = list(pageLength = 9))
  })

  output$regional_heatmap_diff <- renderPlot({
    req(regional_analysis_results())
    regional_analysis_results()$heatmap_diff
  })

  output$regional_heatmap_zscore <- renderPlot({
    req(regional_analysis_results())
    regional_analysis_results()$heatmap_zscore
  })

  output$download_regional_scatter <- downloadHandler(
    filename = function() { "regional_scatter_plot.png" },
    content = function(file) {
      png(file, width = 800, height = 800, res = 120)
      print(regional_analysis_results()$scatter_plot)
      dev.off()
    }
  )

  output$download_regional_heatmap_diff <- downloadHandler(
    filename = function() { "regional_heatmap_diff.png" },
    content = function(file) {
      png(file, width = 1200, height = 1000, res = 120)
      print(regional_analysis_results()$heatmap_diff)
      dev.off()
    }
  )

  output$download_regional_heatmap_zscore <- downloadHandler(
    filename = function() { "regional_heatmap_zscore.png" },
    content = function(file) {
      png(file, width = 1200, height = 1000, res = 120)
      print(regional_analysis_results()$heatmap_zscore)
      dev.off()
    }
  )

  # Venn diagram analysis
  observeEvent(input$run_venn_analysis, {
    req(igp_expr_data(), igp_group_data(), glycan_type_data())
    withProgress(message = 'Running Venn Diagram Analysis...', {
      res <- venn_diagram_analysis(igp_expr_data(), igp_group_data(), glycan_type_data())
      venn_analysis_results(res)
    })
  })

  output$venn_plot_1 <- renderPlot({
    req(venn_analysis_results())
    venn_analysis_results()$venn_plot_1
  })

  output$venn_plot_2 <- renderPlot({
    req(venn_analysis_results())
    venn_analysis_results()$venn_plot_2
  })

  output$venn_plot_3 <- renderPlot({
    req(venn_analysis_results())
    venn_analysis_results()$venn_plot_3
  })

  output$significant_heatmap <- renderPlot({
    req(venn_analysis_results())
    venn_analysis_results()$significant_heatmap
  })

  output$crlm_significant_heatmap <- renderPlot({
    req(venn_analysis_results())
    venn_analysis_results()$crlm_significant_heatmap
  })

  output$crlm_sorted_heatmap <- renderPlot({
    req(venn_analysis_results())
    venn_analysis_results()$crlm_sorted_heatmap
  })

  output$venn_genes_table <- renderDT({
    req(venn_analysis_results())
    genes_df <- data.frame(
      Gene_Type = c("Significant Difference Gene", "CRC vs CRLM Significant Difference Gene", "CRLM vs LIHC Significant Difference Gene", 
                   "Region Method Screening Gene", "Heatmap Method Screening Gene"),
      Count = c(length(venn_analysis_results()$significant_genes),
                length(venn_analysis_results()$final_significant_genes_crc),
                length(venn_analysis_results()$final_significant_genes_lihc),
                length(venn_analysis_results()$different_gene_mean_calculate_regions),
                length(venn_analysis_results()$combined_genes)),
      Genes = c(paste(venn_analysis_results()$significant_genes, collapse = ", "),
                paste(venn_analysis_results()$final_significant_genes_crc, collapse = ", "),
                paste(venn_analysis_results()$final_significant_genes_lihc, collapse = ", "),
                paste(venn_analysis_results()$different_gene_mean_calculate_regions, collapse = ", "),
                paste(venn_analysis_results()$combined_genes, collapse = ", "))
    )
    datatable(genes_df, options = list(scrollX = TRUE, pageLength = 5))
  })

  # Venn diagram download buttons
  output$download_venn_plot_1 <- downloadHandler(
    filename = function() { "venn_plot_1.png" },
    content = function(file) {
      png(file, width = 800, height = 600, res = 120)
      print(venn_analysis_results()$venn_plot_1)
      dev.off()
    }
  )

  output$download_venn_plot_2 <- downloadHandler(
    filename = function() { "venn_plot_2.png" },
    content = function(file) {
      png(file, width = 800, height = 600, res = 120)
      print(venn_analysis_results()$venn_plot_2)
      dev.off()
    }
  )

  output$download_venn_plot_3 <- downloadHandler(
    filename = function() { "venn_plot_3.png" },
    content = function(file) {
      png(file, width = 800, height = 600, res = 120)
      print(venn_analysis_results()$venn_plot_3)
      dev.off()
    }
  )

  output$download_significant_heatmap <- downloadHandler(
    filename = function() { "significant_heatmap.png" },
    content = function(file) {
      png(file, width = 1200, height = 1000, res = 120)
      print(venn_analysis_results()$significant_heatmap)
      dev.off()
    }
  )

  output$download_crlm_significant_heatmap <- downloadHandler(
    filename = function() { "crlm_significant_heatmap.png" },
    content = function(file) {
      png(file, width = 1200, height = 1000, res = 120)
      print(venn_analysis_results()$crlm_significant_heatmap)
      dev.off()
    }
  )

  output$download_crlm_sorted_heatmap <- downloadHandler(
    filename = function() { "crlm_sorted_heatmap.png" },
    content = function(file) {
      png(file, width = 1200, height = 1000, res = 120)
      print(venn_analysis_results()$crlm_sorted_heatmap)
      dev.off()
    }
  )

  output$download_venn_genes <- downloadHandler(
    filename = function() { "venn_genes.csv" },
    content = function(file) {
      genes_df <- data.frame(
        Gene_Type = c("Significant Difference Gene", "CRC vs CRLM Significant Difference Gene", "CRLM vs LIHC Significant Difference Gene", 
                     "Region Method Screening Gene", "Heatmap Method Screening Gene"),
        Count = c(length(venn_analysis_results()$significant_genes),
                  length(venn_analysis_results()$final_significant_genes_crc),
                  length(venn_analysis_results()$final_significant_genes_lihc),
                  length(venn_analysis_results()$different_gene_mean_calculate_regions),
                  length(venn_analysis_results()$combined_genes)),
        Genes = c(paste(venn_analysis_results()$significant_genes, collapse = ", "),
                  paste(venn_analysis_results()$final_significant_genes_crc, collapse = ", "),
                  paste(venn_analysis_results()$final_significant_genes_lihc, collapse = ", "),
                  paste(venn_analysis_results()$different_gene_mean_calculate_regions, collapse = ", "),
                  paste(venn_analysis_results()$combined_genes, collapse = ", "))
      )
      write.csv(genes_df, file, row.names = FALSE)
    }
  )

  # LDA plot and GO analysis
  observeEvent(input$run_lda_go_analysis, {
    req(igp_expr_data(), igp_group_data(), glycan_type_data())
    withProgress(message = 'Running LDA Plot and GO Analysis...', {
      res <- lda_go_analysis(igp_expr_data(), igp_group_data(), glycan_type_data())
      lda_go_analysis_results(res)
    })
  })

  output$pca_plot <- renderPlot({
    req(lda_go_analysis_results())
    lda_go_analysis_results()$pca_plot
  })

  output$lda_plot <- renderPlot({
    req(lda_go_analysis_results())
    lda_go_analysis_results()$lda_plot
  })

  output$svm_plot <- renderPlot({
    req(lda_go_analysis_results())
    lda_go_analysis_results()$svm_plot
  })

  output$go_barplot <- renderPlot({
    req(lda_go_analysis_results())
    if (is.null(lda_go_analysis_results()$go_barplot)) {
      plot.new(); text(0.5, 0.5, "No GO Enrichment Results")
    } else {
      lda_go_analysis_results()$go_barplot
    }
  })

  output$go_dotplot <- renderPlot({
    req(lda_go_analysis_results())
    if (is.null(lda_go_analysis_results()$go_dotplot)) {
      plot.new(); text(0.5, 0.5, "No GO Enrichment Results")
    } else {
      lda_go_analysis_results()$go_dotplot
    }
  })

  output$glycan_pie_plot <- renderPlot({
    req(lda_go_analysis_results())
    lda_go_analysis_results()$glycan_pie_plot
  })

  output$glycan_bar_plot <- renderPlot({
    req(lda_go_analysis_results())
    lda_go_analysis_results()$glycan_bar_plot
  })

  output$go_enrichment_table <- renderDT({
    req(lda_go_analysis_results())
    if (is.null(lda_go_analysis_results()$enrichment_results)) {
      datatable(data.frame(Note = "No GO Enrichment Results"))
    } else {
      datatable(lda_go_analysis_results()$enrichment_results, options = list(scrollX = TRUE, pageLength = 10))
    }
  })

  output$glycan_type_table <- renderDT({
    req(lda_go_analysis_results())
    datatable(lda_go_analysis_results()$glycan_data, options = list(pageLength = 10))
  })

  # LDA plot and GO analysis download buttons
  output$download_pca_plot <- downloadHandler(
    filename = function() { "pca_plot.png" },
    content = function(file) {
      png(file, width = 800, height = 600, res = 120)
      print(lda_go_analysis_results()$pca_plot)
      dev.off()
    }
  )

  output$download_lda_plot <- downloadHandler(
    filename = function() { "lda_plot.png" },
    content = function(file) {
      png(file, width = 800, height = 600, res = 120)
      print(lda_go_analysis_results()$lda_plot)
      dev.off()
    }
  )

  output$download_svm_plot <- downloadHandler(
    filename = function() { "svm_plot.png" },
    content = function(file) {
      png(file, width = 800, height = 600, res = 120)
      print(lda_go_analysis_results()$svm_plot)
      dev.off()
    }
  )

  output$download_go_barplot <- downloadHandler(
    filename = function() { "go_barplot.png" },
    content = function(file) {
      png(file, width = 1000, height = 800, res = 120)
      print(lda_go_analysis_results()$go_barplot)
      dev.off()
    }
  )

  output$download_go_dotplot <- downloadHandler(
    filename = function() { "go_dotplot.png" },
    content = function(file) {
      png(file, width = 1000, height = 800, res = 120)
      print(lda_go_analysis_results()$go_dotplot)
      dev.off()
    }
  )

  output$download_glycan_pie_plot <- downloadHandler(
    filename = function() { "glycan_pie_plot.png" },
    content = function(file) {
      png(file, width = 800, height = 600, res = 120)
      print(lda_go_analysis_results()$glycan_pie_plot)
      dev.off()
    }
  )

  output$download_glycan_bar_plot <- downloadHandler(
    filename = function() { "glycan_bar_plot.png" },
    content = function(file) {
      png(file, width = 800, height = 600, res = 120)
      print(lda_go_analysis_results()$glycan_bar_plot)
      dev.off()
    }
  )

  output$download_go_enrichment <- downloadHandler(
    filename = function() { "go_enrichment_results.csv" },
    content = function(file) {
      write.csv(lda_go_analysis_results()$enrichment_results, file, row.names = FALSE)
    }
  )

  output$download_glycan_type <- downloadHandler(
    filename = function() { "glycan_type_statistics.csv" },
    content = function(file) {
      write.csv(lda_go_analysis_results()$glycan_data, file, row.names = FALSE)
    }
  )

  # KEGG analysis
  observeEvent(input$run_kegg_analysis, {
    req(igp_expr_data(), igp_group_data(), glycan_type_data())
    withProgress(message = 'Running KEGG Analysis...', {
      res <- kegg_analysis(igp_expr_data(), igp_group_data(), glycan_type_data())
      kegg_analysis_results(res)
    })
  })

  output$kegg_barplot <- renderPlot({
    req(kegg_analysis_results())
    if (is.null(kegg_analysis_results()$barplot)) {
      plot.new(); text(0.5, 0.5, "No KEGG Enrichment Results")
    } else {
      kegg_analysis_results()$barplot
    }
  })

  output$kegg_emapplot <- renderPlot({
    req(kegg_analysis_results())
    if (is.null(kegg_analysis_results()$emapplot)) {
      plot.new(); text(0.5, 0.5, "No KEGG Enrichment Results")
    } else {
      kegg_analysis_results()$emapplot
    }
  })

  output$kegg_cnetplot <- renderPlot({
    req(kegg_analysis_results())
    if (is.null(kegg_analysis_results()$cnetplot)) {
      plot.new(); text(0.5, 0.5, "No KEGG Enrichment Results")
    } else {
      kegg_analysis_results()$cnetplot
    }
  })

  output$kegg_table <- renderDT({
    req(kegg_analysis_results())
    if (is.null(kegg_analysis_results()$kegg_table)) {
      datatable(data.frame(Note = "No KEGG Enrichment Results"))
    } else {
      datatable(kegg_analysis_results()$kegg_table, options = list(scrollX = TRUE, pageLength = 10))
    }
  })

  # KEGG combined chart output
  output$kegg_combined_plot <- renderPlot({
    req(kegg_analysis_results())
    if (is.null(kegg_analysis_results()$combined_plot)) {
      plot.new(); text(0.5, 0.5, "No KEGG Enrichment Results")
    } else {
      kegg_analysis_results()$combined_plot
    }
  })

  output$download_kegg_barplot <- downloadHandler(
    filename = function() { "kegg_barplot.png" },
    content = function(file) {
      png(file, width = 1000, height = 800, res = 120)
      print(kegg_analysis_results()$barplot)
      dev.off()
    }
  )
  output$download_kegg_emapplot <- downloadHandler(
    filename = function() { "kegg_emapplot.png" },
    content = function(file) {
      png(file, width = 1000, height = 800, res = 120)
      print(kegg_analysis_results()$emapplot)
      dev.off()
    }
  )
  output$download_kegg_cnetplot <- downloadHandler(
    filename = function() { "kegg_cnetplot.png" },
    content = function(file) {
      png(file, width = 1000, height = 800, res = 120)
      print(kegg_analysis_results()$cnetplot)
      dev.off()
    }
  )
  output$download_kegg_table <- downloadHandler(
    filename = function() { "kegg_enrichment_table.csv" },
    content = function(file) {
      write.csv(kegg_analysis_results()$kegg_table, file, row.names = FALSE)
    }
  )

  # KEGG combined chart download
  output$download_kegg_combined <- downloadHandler(
    filename = function() { "kegg_combined_plot.png" },
    content = function(file) {
      png(file, width = 1200, height = 1000, res = 120)
      print(kegg_analysis_results()$combined_plot)
      dev.off()
    }
  )

  # Volcano plot analysis
  observeEvent(input$run_volcano_analysis, {
    req(igp_expr_data(), igp_group_data(), glycan_type_data())
    withProgress(message = 'Running Volcano Plot Analysis...', {
      res <- volcano_analysis(igp_expr_data(), igp_group_data(), glycan_type_data())
      volcano_analysis_results(res)
    })
  })

  # Volcano plot output
  output$volcano_crc_1 <- renderPlot({
    req(volcano_analysis_results())
    volcano_analysis_results()$volcano_crc_1
  })

  output$volcano_crc_2 <- renderPlot({
    req(volcano_analysis_results())
    volcano_analysis_results()$volcano_crc_2
  })

  output$volcano_lihc_1 <- renderPlot({
    req(volcano_analysis_results())
    volcano_analysis_results()$volcano_lihc_1
  })

  output$volcano_lihc_2 <- renderPlot({
    req(volcano_analysis_results())
    volcano_analysis_results()$volcano_lihc_2
  })

  # Volcano plot combined chart output
  output$volcano_combined_plot <- renderPlot({
    req(volcano_analysis_results())
    volcano_analysis_results()$volcano_combined_plot
  })

  # Protein expression plot
  output$protein_expression_plot <- renderPlot({
    req(volcano_analysis_results(), input$protein_id, input$plot_type)
    protein_key <- paste0(input$protein_id, "_", input$plot_type)
    if (protein_key %in% names(volcano_analysis_results()$protein_plots)) {
      volcano_analysis_results()$protein_plots[[protein_key]]
    }
  })

  # Data table
  output$volcano_data_table <- renderDT({
    req(volcano_analysis_results())
    # Merge two volcano plot data
    crc_data <- volcano_analysis_results()$crlm_vs_crc_data
    crc_data$comparison <- "CRLM vs CRC"
    lihc_data <- volcano_analysis_results()$crlm_vs_lihc_data
    lihc_data$comparison <- "CRLM vs LIHC"
    combined_data <- rbind(crc_data, lihc_data)
    datatable(combined_data, options = list(scrollX = TRUE, pageLength = 15))
  })

  output$cluster_data_table <- renderDT({
    req(volcano_analysis_results())
    datatable(volcano_analysis_results()$cluster_data, options = list(scrollX = TRUE, pageLength = 15))
  })

  # Volcano plot download buttons
  output$download_volcano_crc_1 <- downloadHandler(
    filename = function() { "volcano_crc_1.png" },
    content = function(file) {
      png(file, width = 1000, height = 800, res = 120)
      print(volcano_analysis_results()$volcano_crc_1)
      dev.off()
    }
  )

  output$download_volcano_crc_2 <- downloadHandler(
    filename = function() { "volcano_crc_2.png" },
    content = function(file) {
      png(file, width = 1000, height = 800, res = 120)
      print(volcano_analysis_results()$volcano_crc_2)
      dev.off()
    }
  )

  output$download_volcano_lihc_1 <- downloadHandler(
    filename = function() { "volcano_lihc_1.png" },
    content = function(file) {
      png(file, width = 1000, height = 800, res = 120)
      print(volcano_analysis_results()$volcano_lihc_1)
      dev.off()
    }
  )

  output$download_volcano_lihc_2 <- downloadHandler(
    filename = function() { "volcano_lihc_2.png" },
    content = function(file) {
      png(file, width = 1000, height = 800, res = 120)
      print(volcano_analysis_results()$volcano_lihc_2)
      dev.off()
    }
  )

  output$download_protein_plot <- downloadHandler(
    filename = function() { 
      paste0("protein_", input$protein_id, "_", input$plot_type, ".png") 
    },
    content = function(file) {
      req(volcano_analysis_results(), input$protein_id, input$plot_type)
      protein_key <- paste0(input$protein_id, "_", input$plot_type)
      if (protein_key %in% names(volcano_analysis_results()$protein_plots)) {
        png(file, width = 800, height = 600, res = 120)
        print(volcano_analysis_results()$protein_plots[[protein_key]])
        dev.off()
      }
    }
  )

  output$download_volcano_data <- downloadHandler(
    filename = function() { "volcano_data.csv" },
    content = function(file) {
      req(volcano_analysis_results())
      crc_data <- volcano_analysis_results()$crlm_vs_crc_data
      crc_data$comparison <- "CRLM vs CRC"
      lihc_data <- volcano_analysis_results()$crlm_vs_lihc_data
      lihc_data$comparison <- "CRLM vs LIHC"
      combined_data <- rbind(crc_data, lihc_data)
      write.csv(combined_data, file, row.names = FALSE)
    }
  )

  output$download_cluster_data <- downloadHandler(
    filename = function() { "cluster_data.csv" },
    content = function(file) {
      req(volcano_analysis_results())
      write.csv(volcano_analysis_results()$cluster_data, file, row.names = FALSE)
    }
  )

  # Volcano plot combined chart download
  output$download_volcano_combined <- downloadHandler(
    filename = function() { "volcano_combined_plot.png" },
    content = function(file) {
      png(file, width = 1200, height = 1000, res = 120)
      print(volcano_analysis_results()$volcano_combined_plot)
      dev.off()
    }
  )

  output$venn_combined_plot <- renderPlot({
    req(venn_analysis_results())
    venn_analysis_results()$venn_combined_plot
  })

  output$download_venn_combined <- downloadHandler(
    filename = function() { "venn_combined_plot.png" },
    content = function(file) {
      png(file, width = 1200, height = 1000, res = 120)
      print(venn_analysis_results()$venn_combined_plot)
      dev.off()
    }
  )

  # Data processing reactive values
  processing_result <- reactiveVal(NULL)
  
  observeEvent(input$run_processing, {
    req(input$raw_csv_files)
    file_paths <- input$raw_csv_files$datapath
    file_names <- input$raw_csv_files$name
    # Automatically generate raw.file names, prioritize mixX_repeatY information
    get_rawfile_name <- function(fname) {
      m <- regmatches(fname, regexec("mix([0-9]+)_glyco_4h_3ug_60K_repeat([0-9]+)", fname))
      if (length(m[[1]]) == 3) {
        return(paste0("raw", m[[1]][2], "_", m[[1]][3]))
      } else {
        return(NA)
      }
    }
    raw_file_names <- sapply(file_names, get_rawfile_name)
    if (any(is.na(raw_file_names))) {
      # fallback: raw1, raw2, ...
      raw_file_names <- paste0("raw", seq_along(file_paths))
    }
    withProgress(message = 'Processing data...', {
      res <- process_glycopeptide_files_strict(file_paths, raw_file_names)
      processing_result(res$normalized)
      showNotification("Data processing completed!", type = "message")
    })
  })
  
  output$processing_result_table <- renderDT({
    req(processing_result())
    datatable(head(processing_result(), 100), options = list(scrollX = TRUE, pageLength = 10))
  })
  
  output$download_processing_result <- downloadHandler(
    filename = function() { "all_data_zujian_normal.csv" },
    content = function(file) {
      write.csv(processing_result(), file, row.names = FALSE)
    }
  )
  
  output$processing_hist <- renderPlot({
    req(processing_result())
    ggplot(processing_result(), aes(x = final_scale_intensity)) +
      geom_histogram(bins = 50, fill = "#4F8DFD", color = "black", alpha = 0.7) +
      labs(title = "Normalized Intensity Distribution Histogram", x = "final_scale_intensity", y = "Frequency") +
      theme_classic()
  })
  
  output$processing_boxplot <- renderPlot({
    req(processing_result())
    ggplot(processing_result(), aes(x = raw.file, y = final_scale_intensity, fill = raw.file)) +
      geom_boxplot() +
      labs(title = "Normalized Intensity Boxplot", x = "Raw File (raw.file)", y = "final_scale_intensity") +
      theme_classic() +
      theme(legend.position = "none")
  })
}

# Run application
shinyApp(ui = ui, server = server)