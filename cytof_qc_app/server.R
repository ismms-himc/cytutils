# Required packages:
# install.packages("shiny")
# install.packages("shinydashboard")
# install.packages("shinyFiles")
# install.packages("mclust")
# install.packages("dplyr")
# install.packages("ggplot2")
# install.packages("robustbase")
# install.packages("flowCore")

library(shiny)
library(shinydashboard)
library(mclust)
library(dplyr)
library(himc)
library(flowCore)
library(ggplot2)
library(robustbase)
source("./cytof_toolkit_helper_functions.R")
source("./ui.R")
source_dir("./server_modules")

`%>%` <- dplyr::`%>%`

server <- function(input, output, session) {
  # The default maxRequestSize is 5MB. We override this, arbitrarily choosing 
  # 3GB as the size limit for file uploads.
  options(shiny.maxRequestSize=3000*1024^2) 

  # This variable is used to control the sequence of processes, 
  # ensuring that we only render success and error messages appropriately
  sample_background_control_var <- reactiveValues(
                successful_sample_background_report_completion = FALSE,
                sample_background_report_dir_valid = FALSE,
                sample_background_report_dir_invalid = FALSE,
                is_uploaded_file_type_valid = TRUE,
                fcs_file_import_error = FALSE,
                pre_processing_error = FALSE,
                background_report_generation_error = FALSE,
                aggr_sample_background_report_exists = TRUE,
                aggr_sample_background_report_export_error = FALSE,
                aggr_sample_background_report_export_success = FALSE
                )

  sample_background_file_statuses <- reactiveValues(unsuccessful_fcs_file_import_filenames = "",
                unsuccessful_pre_processing_filenames = "",
                unsuccessful_report_generation_filenames = "",
                successful_sample_background_report_completion_filenames = "",
                sample_background_report_dir = "")

  # The below are sourced from "./server_modules". They handle the logic of 
  # what occurs when a file is uploaded or an actionButton is clicked.
  sample_background_observe_event(input, sample_background_control_var, sample_background_file_statuses)

################### Sample Background App Outputs #############################
  output$sample_background_report_dir_selection_status <- renderUI({
    if (sample_background_control_var$sample_background_report_dir_invalid) {
      p(paste("Error: Invalid target location for sample background report chosen. Please ",
        "try again."),
        style = "color: red; font-size: 14px; margin: 10px;")
    } else if (sample_background_control_var$sample_background_report_dir_valid) {
        p(paste("Target location for sample background report successfully chosen. Proceed by",
         "uploading FCS files."),
        style = "color: green; font-size: 14px; margin: 10px;")
    }
  })

  output$sample_background_invalid_file_type <- renderUI({
    if (!sample_background_control_var$is_uploaded_file_type_valid) {
      p(paste("Error: File(s) of incorrect type were uploaded. Please ",
        "re-upload files, ensuring they are FCS files."),
        style = "color: red; font-size: 14px; margin: 10px;")
    }
  })

  output$sample_background_fcs_file_import_errors <- renderUI({
    error_message <- vector("list")
    if (sample_background_control_var$fcs_file_import_error) {
      for (i in 1:length(sample_background_file_statuses$unsuccessful_fcs_file_import_filenames)) {
        if (i == 1) {
          error_message[[i]] <- list(
            div(style = "color: red; font-size: 14px; margin: 5px; padding-top: 15px;",
              p(paste("Error: The following file(s) failed the file import",
              "step. As such, further processing and data export",
              "did not occur:"))
            )
          )
        } else {
          error_message[[i]] <- list(
            tags$li(style = "color: red; font-size: 14px; margin: 5px; padding-left: 30px;", 
              sample_background_file_statuses$unsuccessful_fcs_file_import_filenames[i])
          )
        }
      }
    }

    error_message
  })

  output$sample_background_pre_processing_errors <- renderUI({
    error_message <- vector("list")
    if (sample_background_control_var$pre_processing_error) {
      for (i in 1:length(sample_background_file_statuses$unsuccessful_pre_processing_filenames)) {
        if (i == 1) {
          error_message[[i]] <- list(
            div(style = "color: red; font-size: 14px; margin: 5px; padding-top: 15px;",
              p(paste("Error: The following files failed the",
              "pre-processing step. As such, further processing",
              "and data export did not occur:"))
            )
          )
        } else {
          error_message[[i]] <- list(
            tags$li(style = "color: red; font-size: 14px; margin: 5px; padding-left: 30px;", 
              sample_background_file_statuses$unsuccessful_pre_processing_filenames[i])
          )
        }
      }
    }

    error_message
  })

  output$sample_background_report_generation_errors <- renderUI({
    error_message <- vector("list")
    if (sample_background_control_var$background_report_generation_error) {
      for (i in 1:length(sample_background_file_statuses$unsuccessful_report_generation_filenames)) {
        if (i == 1) {
          error_message[[i]] <- list(
            div(style = "color: red; font-size: 14px; margin: 5px; padding-top: 15px;",
              p(paste("Error: The following files failed the",
              "background report generating step. As such, data",
              "from these files is not included in the sample background report:"))
            )
          )
        } else {
          error_message[[i]] <- list(
            tags$li(style = "color: red; font-size: 14px; margin: 5px; padding-left: 30px;", 
              sample_background_file_statuses$unsuccessful_report_generation_filenames[i])
          )
        }
      }
    }

    error_message
  })

  output$aggr_sample_background_report_generation_errors <- renderUI({
    if (!sample_background_control_var$aggr_sample_background_report_exists) {
      p(paste("Error: A sample background report was not successfully generated",
        "for any uploaded files and thus, there was no sample background report to export."),
        style = "color: red; font-size: 14px; margin: 10px;")
    }
  })

  output$sample_background_report_export_status <- renderUI({
    if (sample_background_control_var$aggr_sample_background_report_export_error) {
      p(paste("Error: There was an error exporting the sample background report.",
        "Please ensure the chosen report export location exists."),
        style = "color: red; font-size: 14px; margin: 10px;")
    } else if (sample_background_control_var$aggr_sample_background_report_export_success) {
      success_message <- vector("list")

      for (i in 1:length(sample_background_file_statuses$successful_sample_background_report_completion_filenames)) {
          if (i == 1) {
            success_message[[i]] <- list(
              div(style = "color: green; font-size: 14px; margin: 5px; padding-top: 15px;",
                p(paste("Success. Sample background report exported to chosen location.",
                "Report includes data from the following file(s):"))
              )
            )
          } else {
            success_message[[i]] <- list(
              tags$li(style = "color: green; font-size: 14px; margin: 5px; padding-left: 30px;", 
                sample_background_file_statuses$successful_sample_background_report_completion_filenames[i])
            )
          }
      }

      success_message
    }
  })
}