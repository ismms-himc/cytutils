library(flowCore)
library(shinyFiles)
source("./cytof_toolkit_helper_functions.R")

# default gating configuration
event_channel <- "Event_length"
dna_channel <- "191Di"
bead_gates <- dplyr::data_frame(
  channel = c("140Di", "175Di"),
  min = c(5, 5),
  max = c(Inf, Inf)
)
cofactor <- 5

sample_background_observe_event <- function(input, sample_background_control_var, sample_background_file_statuses) {
shinyDirChoose(input, id = "sample_background_report_dir", roots = getVolumes())

  observeEvent(input$sample_background_report_dir, {
    # We reset the reactive values of our sample_background_control_var so that our error 
    # messages fade when the user attempts to re-choose the target qc report directory 
    # and upload FCS files.
    sample_background_control_var$sample_background_report_dir_valid <- FALSE
    sample_background_control_var$sample_background_report_dir_invalid <- FALSE
    sample_background_file_statuses$sample_background_report_dir <- ""

    sample_background_report_dir <- input$sample_background_report_dir
    home <- normalizePath("~")
    sample_background_report_dir_path <- file.path(home, paste(unlist(sample_background_report_dir$path[-1]), collapse = .Platform$file.sep))

    if (file.exists(sample_background_report_dir_path)) {
      sample_background_control_var$sample_background_report_dir_valid <- TRUE
      sample_background_file_statuses$sample_background_report_dir <- sample_background_report_dir_path
    } else {
      sample_background_control_var$sample_background_report_dir_invalid <- TRUE
    }
  })


  observeEvent(input$sample_background_file, {
    if (sample_background_file_statuses$sample_background_report_dir == "") {
      sample_background_control_var$is_output_dir_chosen_before_upload <- FALSE
    } else {
      background_data_file_data_frame <- input$sample_background_file
      num_files_uploaded <- nrow(background_data_file_data_frame)
      if (num_files_uploaded > 0) {

        # We reset the reactive values of our sample_background_control_var so that our error 
        # messages fade when the user attempts to re-upload files.
        sample_background_control_var$is_uploaded_file_type_valid <- TRUE
        sample_background_control_var$successful_sample_background_report_completion <- FALSE
        sample_background_control_var$fcs_file_import_error <- FALSE
        sample_background_control_var$pre_processing_error <- FALSE
        sample_background_control_var$background_report_generation_error <- FALSE
        sample_background_control_var$background_report_generation_success <- FALSE
        sample_background_control_var$sample_background_report_export_success <- FALSE
        sample_background_control_var$sample_background_report_export_error <- FALSE
        sample_background_control_var$is_output_dir_chosen_before_upload <- TRUE

        # We also reset the reactive values of our sample_background_file_statuses variable to
        # ensure error messages include updated file names after a user
        # re-uploads files.
        sample_background_file_statuses$unsuccessful_pre_processing_filenames <- ""
        sample_background_file_statuses$unsuccessful_report_generation_filenames <- ""
        sample_background_file_statuses$unsuccessful_fcs_file_import_filenames <- ""
        sample_background_file_statuses$successful_report_generation_filenames <- ""
        sample_background_file_statuses$successful_report_export_filenames <- ""
        sample_background_file_statuses$unsuccessful_report_export_filenames <- ""
      }

      for (i in 1:num_files_uploaded) {
        background_data_path <- background_data_file_data_frame[i,]$datapath
        # We check to see if uploaded files have FCS file type. We do not 
        # check if "BCKG" is in filename since that is not always the case.
        if (!flowCore::isFCSfile(background_data_path)) {
          sample_background_control_var$is_uploaded_file_type_valid <- FALSE
          return()
        }
      }

      withProgress(value = 0, {

        for (i in 1:num_files_uploaded) {
          incProgress(amount = (1 / num_files_uploaded), 
                message = paste("Processing file", i, "/",
                num_files_uploaded))
          background_data_filename <- background_data_file_data_frame[i,]$name

          background_data_path <- background_data_file_data_frame[i,]$datapath

          background_fcs_data <- fcs_import_file_error_handler(background_data_path)
          if (is.null(background_fcs_data)) {
            sample_background_control_var$fcs_file_import_error <- TRUE
            sample_background_file_statuses$unsuccessful_fcs_file_import_filenames <- c(sample_background_file_statuses$unsuccessful_fcs_file_import_filenames,
                                                                    background_data_filename)
            next
          }

          background_fcs_data_pre_processing <- mass_cytometry_pre_processing_error_handler(background_fcs_data$data, cofactor, 
                              event_channel, dna_channel, 
                              bead_gates)

          if (is.null(background_fcs_data_pre_processing)) {
            sample_background_control_var$pre_processing_error <- TRUE
            sample_background_file_statuses$unsuccessful_pre_processing_filenames <- c(sample_background_file_statuses$unsuccessful_pre_processing_filenames,
                                                                        background_data_filename)

            next
          }

          sample_background_report <- generate_sample_background_report_error_handler(background_fcs_data, 
                                    background_data_filename, 
                                    background_fcs_data_pre_processing,
                                    cofactor)


          if (is.null(sample_background_report)) {
            sample_background_control_var$background_report_generation_error <- TRUE
            sample_background_file_statuses$unsuccessful_report_generation_filenames <- c(sample_background_file_statuses$unsuccessful_report_generation_filenames,
                                                                        background_data_filename)
            next
          } else {
            sample_background_control_var$background_report_generation_success <- TRUE
            sample_background_file_statuses$successful_report_generation_filenames <- c(sample_background_file_statuses$successful_report_generation_filenames, 
                                                                      background_data_filename)
          }

          sample_background_report_export_status <- sample_background_report_export_error_handler(
                                                            sample_background_file_statuses$sample_background_report_dir,
                                                            sample_background_report,
                                                            background_data_filename)

          if (is.null(sample_background_report_export_status)) {
            sample_background_control_var$sample_background_report_export_success <- TRUE
            sample_background_file_statuses$successful_report_export_filenames <- c(sample_background_file_statuses$successful_report_export_filenames,
                                                                                    background_data_filename)
            sample_background_control_var$sample_background_report_tables[[background_data_filename]] <- sample_background_report
          } else if (sample_background_report_export_status == "Unsuccessful"){
            sample_background_control_var$sample_background_report_export_error <- TRUE
            sample_background_file_statuses$unsuccessful_report_export_filenames <- c(sample_background_file_statuses$unsuccessful_report_export_filenames,
                                                                                      background_data_filename)
          }          

          rm(background_fcs_data_pre_processing,
              sample_background_report)
        }
      })
    }
  })
}