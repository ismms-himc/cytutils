library(shiny)
library(shinydashboard)
library(shinyFiles)

qc_report_app_tab_item <- tabItem(tabName = "qc_report_app",
                      h2("CyTOF QC Report Generator",
                        style = "padding-left: 15px"),
                      fluidRow(
                        sidebarPanel(
                          p(paste("Please select the location you would like to", 
                          "export the QC report to. A message will appear to the",
                          "right when the location has successfully been",
                          "selected.")),
                          shinyDirButton(id="cytof_qc_report_dir", 
                                        label="Choose directory", 
                                        title="Select directory" ),
                          br(),
                          br(),
                          br(),
                          p(paste("Please upload CyTOF FCS file(s) to run a QC report.",
                                    "Report results will be exported in a",
                                    "cytof_qc_report.csv in the location chosen",
                                    "above.")),
                          br(),
                          checkboxGroupInput(inputId = "checkGroup", 
                                             label = NULL, 
                                             choices = list("Generate Channel vs Time plot" = 1)
                                             ),
                          fileInput(
                            # NB: Regarding all inputId's, inputId is assigned to 
                            # singular "x_file" because under most circumstances, a single
                            # file will be uploaded; however, multiple files may be 
                            # uploaded. Whenever a file upload completes, the corresponding 
                            # input variable is set to a dataframe. This dataframe contains 
                            # one row for each selected file. 
                            inputId = "cytof_qc_file",
                            label = paste("Please do not exit the application until a status", 
                                          "message appears to the right of this",
                                          "section or data export may be compromised."),
                            multiple = TRUE,
                            accept = c(".fcs", "fcs")
                          )
                        ),
                        mainPanel(
                          htmlOutput("cytof_qc_report_dir_selection_status"),
                          htmlOutput("cytof_qc_invalid_file_type"),
                          htmlOutput("cytof_qc_fcs_file_import_errors"),
                          htmlOutput("cytof_qc_pre_processing_errors"),
                          htmlOutput("cytof_qc_report_generation_errors"),
                          htmlOutput("cytof_qc_export_successes"),
                          htmlOutput("cytof_qc_export_errors")
                        )
                      ),
                      uiOutput("gating_inspection_and_visualization")
                  )