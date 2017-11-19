# Configuration:
# install.packages("shiny")
# install.packages("shinydashboard")
# install.packages("shinyFiles")
# install_github("ismmshimc/cytutils")

library(shiny)
library(shinydashboard)
library(shinyFiles)
library(cytutils)
source("./helper_functions.R")
source("./ui.R")
source_dir("./server_modules")

`%>%` <- dplyr::`%>%`

server <- function(input, output, session) {
  # The default maxRequestSize is 5MB. We override this, arbitrarily choosing 
  # 3GB as the size limit for file uploads.
  options(shiny.maxRequestSize=3000*1024^2) 

  # The "x_control_var" variables are used to control the sequence of processes, 
  # ensuring that we only render a message that the QC report data was uploaded 
  # to Google Drive after we use the googlesheets R package to do so.
  fcs_channel_rename_control_var <- reactiveValues(successful_completion = FALSE,
                is_fcs_file_dir_chosen = FALSE,
                fcs_file_detection_error = FALSE
              )
  fcs_channel_rename_statuses <- reactiveValues(
                fcs_channel_rename_dir_path = ""
              )



  # NB: The below are sourced from "./server_modules". They handle the logic of 
  # what occurs when a file is uploaded or an actionButton is clicked.
  fcs_channel_rename_observe_event(input, fcs_channel_rename_control_var, fcs_channel_rename_statuses)

  ################### CyTOF QC Report App Outputs #############################
  # This occurs when the user chooses a directory with FCS files to rename.
  output$new_channel_rename_csv_generation_status <- renderUI({
    if (fcs_channel_rename_control_var$is_fcs_file_dir_chosen) {
      # TODO: run cytutils::channel_rename with fcs_channel_rename_statuses$fcs_channel_rename_dir_path
      # as an argument
      # p(paste("Success: File(s) of incorrect type were uploaded. This ",
      #   "application only supports FCS files."),
      #   style = "color: red; font-size: 14px; margin: 10px;")
    } else if (fcs_channel_rename_control_var$fcs_file_detection_error) {
      p("Error: FCS files were not detected in the selected directory."),
        style = "color: red; font-size: 14px; margin: 10px;")
    }
  })

  # output$import_channel_rename_csv_and_export_fcs_files_status <- renderUI({
  #   if () {
  #     p(paste("Error:...................", 
  #       ".............."),
  #       style = "color: red; font-size: 14px; margin: 10px;")
  #   }
  # })
}