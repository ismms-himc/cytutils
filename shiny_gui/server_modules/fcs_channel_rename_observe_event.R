library(tools)

fcs_channel_rename_observe_event <- function(input, fcs_file_dir_control_var, fcs_file_dir_statuses, fcs_file_selection_control_var) {
  shinyDirChoose(input, id = "fcs_channel_rename_dir", roots = c(home = '~'), filetypes = c('', 'fcs'))
  

  observeEvent(input$fcs_channel_rename_dir, {

    fcs_channel_rename_dir <- input$fcs_channel_rename_dir
    home <- normalizePath("~")
    fcs_channel_rename_dir_path <- file.path(home, paste(unlist(fcs_channel_rename_dir$path[-1]), collapse = .Platform$file.sep))
    files <- list.files(fcs_channel_rename_dir)

    # Storing FCS file names is not necessary to support current functionality.
    # We store them in the event that we wish to render more specific success/error
    # messages in future app versions.
    fcs_files_detected <- c()

    for (file in files) {
      if (file_ext(file) == "fcs") {
        fcs_files_detected <- c(fcs_files_detected, file)
      }
    }

    if (length(fcs_files_detected > 0)) {
      # TODO: add status messages
      # We update/reset the reactive values of our control_var and statuses so that our old status 
      # messages fade when the user attempts to re-select a directory.
      fcs_channel_rename_statuses$fcs_channel_rename_dir_path <- fcs_channel_rename_dir_path

      fcs_channel_rename_control_var$is_fcs_file_dir_chosen <- TRUE
      fcs_channel_rename_control_var$fcs_file_detection_error <- FALSE
    } else {
      fcs_channel_rename_statuses$fcs_channel_rename_dir_path <- ""

      fcs_channel_rename_control_var$is_fcs_file_dir_chosen <- FALSE
      fcs_channel_rename_control_var$fcs_file_detection_error <- TRUE
    }
  })
}