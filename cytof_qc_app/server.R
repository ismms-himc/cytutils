# Required packages:
# install.packages("shiny")
# install.packages("shinydashboard")
# install.packages("shinyFiles")
# install.packages("mclust")
# install.packages("dplyr")
# install.packages("ggplot2")
# install.packages("robustbase")
# install.packages("flowCore")
# install.packages("matrixStats")

library(shiny)
library(shinydashboard)
library(mclust)
library(dplyr)
library(flowCore)
library(ggplot2)
library(robustbase)
library(matrixStats)
source("./cytof_toolkit_helper_functions.R")
source("./ui.R")
source_dir("./server_modules")

`%>%` <- dplyr::`%>%`

server <- function(input, output, session) {
  # The default maxRequestSize is 5MB. We override this, arbitrarily choosing 
  # 3GB as the size limit for file uploads.
  options(shiny.maxRequestSize=3000*1024^2) 

  cytof_qc_control_var <- reactiveValues(cytof_qc_report_dir_valid = FALSE,
                cytof_qc_report_dir_invalid = FALSE,
                is_output_dir_chosen_before_upload = TRUE,
                # aggr_cytof_qc_report_export_error = FALSE,
                successful_completion = FALSE,
                is_uploaded_file_type_valid = TRUE,
                fcs_file_import_error = FALSE,
                pre_processing_error = FALSE,
                qc_report_generation_error = FALSE,
                render_gating_inspection = FALSE,
                abnormal_gating_flag_error = FALSE,
                successful_abnormal_gating_flag = FALSE,
                abnormal_gating_unflag_error = FALSE,
                successful_abnormal_gating_unflag = FALSE,
                render_manual_gating = FALSE,
                manual_gating_error = FALSE,
                manual_gating_success = FALSE,
                updated_qc_report_generation_error = FALSE
              )

  cytof_qc_file_statuses <- reactiveValues(cytof_qc_report_dir = "",
                unsuccessful_fcs_file_import_filenames = "",
                unsuccessful_pre_processing_filenames = "",
                unsuccessful_report_generation_filenames = "",
                successful_completion_filenames = "",
                successful_abnormal_gating_flag_filename = "",
                unsuccessful_abnormal_gating_flag_filename = "",
                unsuccessful_abnormal_gating_unflag_filename = "",
                successful_abnormal_gating_unflag_filename = "",
                unsuccessful_updated_qc_report_filename = ""
              )

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
                background_report_generation_success = FALSE,
                sample_background_report_export_success = FALSE,
                sample_background_report_export_error = FALSE,
                is_output_dir_chosen_before_upload = TRUE
                )

  sample_background_file_statuses <- reactiveValues(unsuccessful_fcs_file_import_filenames = "",
                unsuccessful_pre_processing_filenames = "",
                unsuccessful_report_generation_filenames = "",
                successful_report_generation_filenames = "",
                sample_background_report_dir = "",
                successful_report_export_filenames = "",
                unsuccessful_report_export_filenames = "")

  # NB: We initialize a list, arbitrarily choosing a length of 1000 (representing
  # a unrealistically high number of files for which qc report data was 
  # successfully exported) in order to allow us to manipulate the 
  # list in amortized constant time.
  # https://stackoverflow.com/questions/2436688/append-an-object-to-a-list-in-r-in-amortized-constant-time-o1
  cytof_qc_gating_inspection <- reactiveValues(pre_processed_data = vector("list", 1000),
                                              currently_rendered_gating_filename = "")

  # The below are sourced from "./server_modules". They handle the logic of 
  # what occurs when a file is uploaded or an actionButton is clicked.
  cytof_qc_observe_event(input, cytof_qc_control_var, cytof_qc_file_statuses, cytof_qc_gating_inspection)
  sample_background_observe_event(input, sample_background_control_var, sample_background_file_statuses)
################### CyTOF QC Feature Outputs #############################

  output$cytof_qc_report_dir_selection_status <- renderUI({
    if (!cytof_qc_control_var$is_output_dir_chosen_before_upload) {
        p(paste("Target location for cytof qc report must be chosen prior",
         "to uploading FCS files."),
        style = "color: red; font-size: 14px; margin: 10px;")
    } else if (cytof_qc_control_var$cytof_qc_report_dir_invalid) {
      p(paste("Error: Invalid target location for cytof qc report chosen. Please ",
        "try again."),
        style = "color: red; font-size: 14px; margin: 10px;")
    } else if (cytof_qc_control_var$cytof_qc_report_dir_valid) {
        p(paste("Target location for cytof qc report successfully chosen. Proceed by",
         "uploading FCS files."),
        style = "color: green; font-size: 14px; margin: 10px;")
    } 
  })

































































































































################### Sample Background Feature Outputs #############################
  output$sample_background_report_dir_selection_status <- renderUI({
    if (!sample_background_control_var$is_output_dir_chosen_before_upload) {
        p(paste("Target location for sample background report must be chosen prior",
         "to uploading FCS files."),
        style = "color: red; font-size: 14px; margin: 10px;")
    } else if (sample_background_control_var$sample_background_report_dir_invalid) {
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

  output$sample_background_report_export_errors <- renderUI({
    error_message <- vector("list")
    if (sample_background_control_var$sample_background_report_export_error) {
      for (i in 1:length(sample_background_file_statuses$unsuccessful_report_export_filenames)) {
        if (i == 1) {
          error_message[[i]] <- list(
            div(style = "color: red; font-size: 14px; margin: 5px; padding-top: 15px;",
              p(paste("Error: Sample background report export failed for the",
              "following files:"))
            )
          )
        } else {
          error_message[[i]] <- list(
            tags$li(style = "color: red; font-size: 14px; margin: 5px; padding-left: 30px;", 
              sample_background_file_statuses$unsuccessful_report_export_filenames[i])
          )
        }
      }
    }

    error_message
  })

  output$sample_background_report_export_successes <- renderUI({
    success_message <- vector("list")
    if (sample_background_control_var$sample_background_report_export_success) {
      for (i in 1:length(sample_background_file_statuses$successful_report_export_filenames)) {
        if (i == 1) {
          success_message[[i]] <- list(
            div(style = "color: green; font-size: 14px; margin: 5px; padding-top: 15px;",
              p(paste("Success: Sample background reports were exported for the",
              "following files:"))
            )
          )
        } else {
          success_message[[i]] <- list(
            tags$li(style = "color: green; font-size: 14px; margin: 5px; padding-left: 30px;", 
              sample_background_file_statuses$successful_report_export_filenames[i])
          )
        }
      }
    }

    success_message
  })
}