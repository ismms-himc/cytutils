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
                updated_qc_report_generation_error = FALSE,
                qc_report_export_success = FALSE,
                qc_report_export_error = FALSE
              )

  cytof_qc_file_statuses <- reactiveValues(cytof_qc_report_dir = "",
                unsuccessful_fcs_file_import_filenames = "",
                unsuccessful_pre_processing_filenames = "",
                unsuccessful_report_generation_filenames = "",
                # successful_completion_filenames = "",
                successful_abnormal_gating_flag_filename = "",
                unsuccessful_abnormal_gating_flag_filename = "",
                unsuccessful_abnormal_gating_unflag_filename = "",
                successful_abnormal_gating_unflag_filename = "",
                unsuccessful_updated_qc_report_filename = "",
                successful_report_export_filenames = "",
                unsuccessful_report_export_filenames = ""
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
  cytof_qc_gating_visualization_observe_event(input, output, session, cytof_qc_control_var, cytof_qc_gating_inspection)
  flag_abnormal_gating_observe_event(input, cytof_qc_gating_inspection, cytof_qc_control_var, cytof_qc_file_statuses)
  # unflag_abnormal_gating_observe_event(input, cytof_qc_control_var, cytof_qc_file_statuses)
  # cytof_qc_manually_update_gating_observe_event(input, output, session, cytof_qc_control_var, cytof_qc_file_statuses, cytof_qc_gating_inspection)
  # cytof_qc_generate_updated_qc_report_observe_event(input, output, cytof_qc_control_var, cytof_qc_file_statuses, cytof_qc_gating_inspection)
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

  output$cytof_qc_invalid_file_type <- renderUI({
    if (!cytof_qc_control_var$is_uploaded_file_type_valid) {
      p(paste("Error: File(s) of incorrect type were uploaded. This ",
        "application only supports FCS files."),
        style = "color: red; font-size: 14px; margin: 10px;")
    }
  })


  output$cytof_qc_fcs_file_import_errors <- renderUI({
    error_message <- vector("list")
    if (cytof_qc_control_var$fcs_file_import_error) {
      for (i in 1:length(cytof_qc_file_statuses$unsuccessful_fcs_file_import_filenames)) {
        if (i == 1) {
          error_message[[i]] <- list(
            div(style = "color: red; font-size: 14px; margin: 5px; padding-top: 15px;",
              p(paste("Error: The following file(s) failed the file import",
              "step. As such, pre-processing, QC report generation, and data",
              "transfer to Google Drive did not occur:"))
            )
          )
        } else {
          error_message[[i]] <- list(
            tags$li(style = "color: red; font-size: 14px; margin: 5px; padding-left: 30px;", 
              cytof_qc_file_statuses$unsuccessful_fcs_file_import_filenames[i])
          )
        }
      }
    }

    error_message
  })

  output$cytof_qc_pre_processing_errors <- renderUI({
    error_message <- vector("list")
    if (cytof_qc_control_var$pre_processing_error) {
      for (i in 1:length(cytof_qc_file_statuses$unsuccessful_pre_processing_filenames)) {
        if (i == 1) {
          error_message[[i]] <- list(
            div(style = "color: red; font-size: 14px; margin: 5px; padding-top: 15px;",
              p(paste("Error: The following file(s) failed the",
              "pre-processing step. As such, QC report generation and data",
              "export did not occur:"))
            )
          )
        } else {
          error_message[[i]] <- list(
            tags$li(style = "color: red; font-size: 14px; margin: 5px; padding-left: 30px;", 
              cytof_qc_file_statuses$unsuccessful_pre_processing_filenames[i])
          )
        }
      }
    }

    error_message
  })

  output$cytof_qc_report_generation_errors <- renderUI({
    error_message <- vector("list")
    if (cytof_qc_control_var$qc_report_generation_error) {
      for (i in 1:length(cytof_qc_file_statuses$unsuccessful_report_generation_filenames)) {
        if (i == 1) {
          error_message[[i]] <- list(
            div(style = "color: red; font-size: 14px; margin: 5px; padding-top: 15px;",
              p(paste("Error: The following file(s) failed the",
              "QC report generation step. As such, data export",
              "did not occur:"))
            )
          )
        } else {
          error_message[[i]] <- list(
            tags$li(style = "color: red; font-size: 14px; margin: 5px; padding-left: 30px;", 
              cytof_qc_file_statuses$unsuccessful_report_generation_filenames[i])
          )
        }
      }
    }

    error_message
  })

  output$cytof_qc_export_errors <- renderUI({
    error_message <- vector("list")
    if (cytof_qc_control_var$qc_report_export_error) {
      for (i in 1:length(cytof_qc_file_statuses$unsuccessful_report_export_filenames)) {
        if (i == 1) {
          error_message[[i]] <- list(
            div(style = "color: red; font-size: 14px; margin: 5px; padding-top: 15px;",
              p(paste("Error: Data from the following file(s)",
              "was not exported. Please ensure target folder exists."))
            )
          )
        } else {
          error_message[[i]] <- list(
            tags$li(style = "color: red; font-size: 14px; margin: 5px; padding-left: 30px;", 
              cytof_qc_file_statuses$unsuccessful_report_export_filenames[i])
          )
        }
      }
    }

    error_message
  })

  output$cytof_qc_export_successes <- renderUI({
    success_message <- vector("list")
    if (cytof_qc_control_var$qc_report_export_success) {
      for (i in 1:length(cytof_qc_file_statuses$successful_report_export_filenames)) {
        if (i == 1) {
          success_message[[i]] <- list(
            div(style = "color: green; font-size: 14px; margin: 5px; padding-top: 15px;",
              p(paste("Success. QC report data from the uploaded",
              "file(s) has been exported to the selected output folder:"))
            )
          )
        } else {
          success_message[[i]] <- list(
            tags$li(style = "color: green; font-size: 14px; margin: 5px; padding-left: 30px;", 
              cytof_qc_file_statuses$successful_report_export_filenames[i])
          )
        }
      }
    }

    success_message
  })

  output$gating_inspection_and_visualization <- renderUI({
    if (cytof_qc_control_var$render_gating_inspection) {
      fluidRow(
        tags$hr(style="border-color: #C0C0C0;"),
        column(4,
          h3("Gating Inspection"),
          p("Please choose a file name to render a gating visualization:"),
          selectInput(inputId = "gating_filename", 
                      label = "Filename",
                      choices = cytof_qc_file_statuses$successful_report_export_filenames),
          actionButton(inputId = "generate_gating_visualization",
                      label = "Generate Gating Visualization", 
                      icon = icon("area-chart")),
          uiOutput("abnormal_gating_flag"),
          htmlOutput("abnormal_gating_flag_status_message"),
          uiOutput("undo_abnormal_gating_flag"),
          htmlOutput("undo_abnormal_gating_flag_status_message")
        ),
        column(5,
          plotOutput("gating_visualization",
            brush = brushOpts(
              id = "manual_gating_brush"
              ))
        ), 
        column(3,
          uiOutput("manual_gating"),
          htmlOutput("update_gating_status_message"),
          uiOutput("update_qc_report"),
          htmlOutput("generate_updated_qc_report_error_message")
          # htmlOutput("transfer_updated_qc_report_to_google_drive_status_message")
        )
      )     
    }
  })

  output$abnormal_gating_flag_status_message <- renderUI({
    status_message <- vector("list")
    if (cytof_qc_control_var$successful_abnormal_gating_flag & !cytof_qc_control_var$abnormal_gating_unflag_error) {
      status_message[[1]] <- list(
            div(style = "color: green; font-size: 12px; margin: 5px; padding-top: 15px;",
              p(paste("Success. The following file was flagged as having abnormal",
              "gating in the exported QC report:")
              )
            )
          )

      status_message[[2]] <- list(
          tags$li(style = "color: green; font-size: 12px; margin: 5px; padding-left: 30px;", 
              cytof_qc_file_statuses$successful_abnormal_gating_flag_filename)
        )
    } else if (cytof_qc_control_var$abnormal_gating_flag_error) {
      status_message[[1]] <- list(
            div(style = "color: red; font-size: 12px; margin: 5px; padding-top: 15px;",
              p(paste("Error: The following file was not flagged as having abnormal",
              "gating in the exported QC report. Please ensure the previously chosen",
              "target directory still exists and the previously exported QC report's",
              "filename was not changed.")
              )
            )
          )

      status_message[[2]] <- list(
          tags$li(style = "color: red; font-size: 12px; margin: 5px; padding-left: 30px;", 
              cytof_qc_file_statuses$unsuccessful_abnormal_gating_flag_filename)
        )
    }

    status_message
  })

  # # We render an undo button when abnormal gating was successfully flagged.
  # output$undo_abnormal_gating_flag <- renderUI({
  #   if (cytof_qc_control_var$successful_abnormal_gating_flag) {
  #     actionButton(inputId = "unflag_abnormal_gating",
  #         label = "Undo Abnormal Gating Flag", 
  #         icon = icon("undo"))
  #   }
  # })

  # output$undo_abnormal_gating_flag_status_message <- renderUI({
  #   status_message <- vector("list")
  #   if (cytof_qc_control_var$successful_abnormal_gating_unflag) {
  #     status_message[[1]] <- list(
  #           div(style = "color: green; font-size: 12px; margin: 5px; padding-top: 15px;",
  #             p(paste("Success. The abnormal gating flag for the following file",
  #             "was removed from Google Drive:")
  #             )
  #           )
  #         )

  #     status_message[[2]] <- list(
  #         tags$li(style = "color: green; font-size: 12px; margin: 5px; padding-left: 30px;", 
  #             cytof_qc_file_statuses$successful_abnormal_gating_unflag_filename)
  #       )
  #   } else if (cytof_qc_control_var$abnormal_gating_unflag_error) {
  #     status_message[[1]] <- list(
  #           div(style = "color: red; font-size: 12px; margin: 5px; padding-top: 15px;",
  #             p(paste("Error: The abnormal gating flag for the following file was not",
  #             "removed from Google Drive. Please ensure the \"acq file path\" column",
  #             "cell corresponding to this file is properly filled in. A modified",
  #             "or missing column name may have also caused the error.")
  #             )
  #           )
  #         )

  #     status_message[[2]] <- list(
  #         tags$li(style = "color: red; font-size: 12px; margin: 5px; padding-left: 30px;", 
  #             cytof_qc_file_statuses$unsuccessful_abnormal_gating_unflag_filename)
  #       )
  #   }

  #   status_message
  # })

  # # We a section to allow the user to manually adjust gating after the user
  # # successfully flags gating as abnormal.
  # output$manual_gating <- renderUI({
  #   if (cytof_qc_control_var$render_manual_gating) {
  #     div(
  #       h4("Manually Adjust Abnormal Gating"),
  #       p(paste("Click and drag on the plot to gate a population. Select a",
  #       "category to assign the population to:")),
  #       selectInput(inputId = "manual_gating_category",
  #                   label = "Reclassify Selection As:",
  #                   choices = c("bead", "cell", "debris")),
  #       actionButton(inputId = "manually_update_gating",
  #                   label = "Update Gating",
  #                   icon = icon("wrench"))  
  #     )
  #   }
  # })

  # output$update_gating_status_message <- renderUI({
  #   # We check if cytof_qc_control_var$render_manual_gating to avoid rendering status
  #   # messages in the case that a user generated a gating visualization for a different
  #   # file and did not yet flag it as having abnormal gating.
  #   if (cytof_qc_control_var$render_manual_gating){
  #     if (cytof_qc_control_var$manual_gating_error) {
  #       p(paste("Error: Gating not updated. Please click and drag on the plot to",
  #         "gate a population. Please also ensure the desired population is",
  #         "completely within manually created gate boundaries."),
  #         style = "color: red; font-size: 12px; margin: 5px; padding-top: 15px;"
  #       )

  #     } else if (cytof_qc_control_var$manual_gating_success) {
  #         p("Success: Gating updated.", style = "color: green; font-size: 12px; margin: 5px; padding-top: 15px;")
  #     }
  #   }
  # })

  # output$update_qc_report <- renderUI({
  #   # We check if cytof_qc_control_var$render_manual_gating to avoid rendering
  #   # a prompt for an updated QC report in case a user generates a gating 
  #   # visualization for a different file without flagging it as having abnormal gating.
  #   if (cytof_qc_control_var$render_manual_gating & cytof_qc_control_var$manual_gating_success) {
  #     div(
  #       tags$hr(style = "border-color: #C0C0C0; margin: 35px 0px 27px 0px"),
  #       h4("Satisfactory Gating?"),
  #         p("If yes, please click below to generate an updated QC report and send",
  #           "it to Google Drive. Otherwise, please continue to update gating."),
  #         actionButton(inputId = "generate_updated_qc_report",
  #               label = "Update QC Report", 
  #               icon = icon("check-circle"),
  #               style = "margin-top: 7px;")
  #     )
  #   }
  # })

  # output$generate_updated_qc_report_error_message <- renderUI({
  #   status_message <- vector("list")

  #   if (cytof_qc_control_var$render_manual_gating & cytof_qc_control_var$updated_qc_report_generation_error) {
  #     status_message[[1]] <- list(
  #           div(style = "color: red; font-size: 12px; margin: 5px; padding-top: 15px;",
  #             p(paste("Error: An updated QC report was not generated for the",
  #             "following file:")
  #             )
  #           )
  #         )
  #     status_message[[2]] <- list(
  #         tags$li(style = "color: red; font-size: 12px; margin: 5px; padding-left: 30px;", 
  #             cytof_qc_file_statuses$unsuccessful_updated_qc_report_filename)
  #       )
  #   }

  #   status_message
  # })

  # output$transfer_updated_qc_report_to_google_drive_status_message <- renderUI({
  #   status_message <- vector("list")
  #   # We check if cytof_qc_control_var$render_manual_gating to avoid rendering status
  #   # messages in the case that a user generated a gating visualization for a different
  #   # file and did not yet flag it as having abnormal gating.
  #   if (cytof_qc_control_var$render_manual_gating){
  #     if (cytof_qc_control_var$google_drive_updated_qc_report_transfer_error) {
  #       status_message[[1]] <- list(
  #           div(style = "color: red; font-size: 12px; margin: 5px; padding-top: 15px;",
  #             p(paste("Error: An updated QC report was not sent to Google Drive",
  #             "for the following file:")
  #             )
  #           )
  #         )
  #     status_message[[2]] <- list(
  #         tags$li(style = "color: red; font-size: 12px; margin: 5px; padding-left: 30px;", 
  #             cytof_qc_file_statuses$unsuccessful_google_drive_updated_qc_report_transfer_filename)
  #       )
  #     } else if (cytof_qc_control_var$google_drive_updated_qc_report_transfer_success) {
  #        status_message[[1]] <- list(
  #           div(style = "color: green; font-size: 12px; margin: 5px; padding-top: 15px;",
  #             p(paste("Success. An updated QC report was generated and sent to",
  #             "Google Drive for the following file:")
  #             )
  #           )
  #         )

  #       status_message[[2]] <- list(
  #           tags$li(style = "color: green; font-size: 12px; margin: 5px; padding-left: 30px;", 
  #               cytof_qc_file_statuses$successful_google_drive_updated_qc_report_transfer_filename)
  #         )
  #     }
  #   }

  #   status_message
  # })































































































































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