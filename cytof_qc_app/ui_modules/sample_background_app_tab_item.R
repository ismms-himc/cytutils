library(shiny)
library(shinydashboard)
library(shinyFiles)

sample_background_app_tab_item <-  tabItem(tabName = "sample_background_app",
      h2("Sample Background Tracking",
        style = "padding-left: 15px"),
      sidebarPanel(
          p(paste("Please select the location you would like to export the sample",
                    "background report to. A message will appear to the right when",
                    "the location has successfully been selected.")),
          shinyDirButton(id="sample_background_report_dir", label="Choose directory", title="Select directory" ),
          br(),
          br(),
          br(),
          p(paste(
            "Please upload sample background FCS file(s). Data from the",
            "uploaded file(s) will be pre-processed and exported.")),
          fileInput(
            inputId = "sample_background_file",
            label = paste("Please do not exit the application until a sample background report-specific status",
              "message appears to the right of this section or QC report export",
              "may be compromised."),
            multiple = TRUE,
            accept = c(".fcs", "fcs")
          )
        ),
        mainPanel(
          htmlOutput("sample_background_report_dir_selection_status"),
          htmlOutput("sample_background_invalid_file_type"),
          htmlOutput("sample_background_fcs_file_import_errors"),
          htmlOutput("sample_background_pre_processing_errors"),
          htmlOutput("sample_background_report_generation_errors"),
          htmlOutput("sample_background_report_successful_completions"),
          htmlOutput("sample_background_report_export_successes"),
          htmlOutput("sample_background_report_export_errors")
        )
      )