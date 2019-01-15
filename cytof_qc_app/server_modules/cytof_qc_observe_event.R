library(flowCore)
library(shinyFiles)

# default gating configuration
event_channel <- "Event_length"
dna_channel <- "191Di"
bead_gates <- dplyr::data_frame(
  channel = c("140Di", "175Di"),
  min = c(5, 5),
  max = c(Inf, Inf)
)
cofactor <- 5

cytof_qc_observe_event <- function(input, cytof_qc_control_var, cytof_qc_file_statuses, cytof_qc_gating_inspection) {
  shinyDirChoose(input, id = "cytof_qc_report_dir", roots = getVolumes())

  observeEvent(input$cytof_qc_report_dir, {
    # We reset the reactive values of our cytof_qc_control_var so that our error 
    # messages fade when the user attempts to re-choose the target qc report directory 
    # and upload FCS files.
    cytof_qc_control_var$cytof_qc_report_dir_valid <- FALSE
    cytof_qc_control_var$cytof_qc_report_dir_invalid <- FALSE
    cytof_qc_file_statuses$cytof_qc_report_dir <- ""

    cytof_qc_report_dir <- input$cytof_qc_report_dir

    if (class(cytof_qc_report_dir) != "integer") {
      cytof_qc_report_dir_path <- file.path(paste(unlist(cytof_qc_report_dir$path[-1]), collapse = .Platform$file.sep))

      
      cytof_qc_report_dir_path <- paste0(.Platform$file.sep, cytof_qc_report_dir_path)

      if (file.exists(cytof_qc_report_dir_path)) {
        cytof_qc_control_var$cytof_qc_report_dir_valid <- TRUE
        cytof_qc_file_statuses$cytof_qc_report_dir <- cytof_qc_report_dir_path
      } else {
        cytof_qc_control_var$cytof_qc_report_dir_invalid <- TRUE
      }
    }
  })

  observeEvent(input$cytof_qc_file, {
    if (cytof_qc_file_statuses$cytof_qc_report_dir == "") {
      cytof_qc_control_var$is_output_dir_chosen_before_upload <- FALSE
    } else {
      fcs_file_data_frame <- input$cytof_qc_file
      num_files_uploaded <- nrow(fcs_file_data_frame)

      if (num_files_uploaded > 0) {
        # We reset the reactive values of our cytof_qc_control_var so that our error 
        # messages fade when the user attempts to re-upload files.
        cytof_qc_control_var$is_uploaded_file_type_valid <- TRUE
        cytof_qc_control_var$fcs_file_import_error <- FALSE
        cytof_qc_control_var$pre_processing_error <- FALSE
        cytof_qc_control_var$qc_report_generation_error <- FALSE
        cytof_qc_control_var$render_gating_inspection <- FALSE
        cytof_qc_control_var$abnormal_gating_flag_error <- FALSE
        cytof_qc_control_var$successful_abnormal_gating_flag <- FALSE
        cytof_qc_control_var$abnormal_gating_unflag_error <- FALSE
        cytof_qc_control_var$successful_abnormal_gating_unflag <- FALSE
        cytof_qc_control_var$render_manual_gating <- FALSE
        cytof_qc_control_var$manual_gating_success <- FALSE
        cytof_qc_control_var$manual_gating_error <- FALSE
        cytof_qc_control_var$updated_qc_report_generation_error <- FALSE
        cytof_qc_control_var$qc_report_export_success <- FALSE
        cytof_qc_control_var$qc_report_export_error <- FALSE
        
        # We also reset the reactive values of our cytof_qc_file_statuses variable to
        # ensure error messages include updated file names after a user
        # re-uploads files.
        cytof_qc_file_statuses$unsuccessful_fcs_file_import_filenames <- ""
        cytof_qc_file_statuses$unsuccessful_pre_processing_filenames <- ""
        cytof_qc_file_statuses$unsuccessful_report_generation_filenames <- ""
        cytof_qc_file_statuses$successful_abnormal_gating_flag_filename <- ""
        cytof_qc_file_statuses$unsuccessful_abnormal_gating_flag_filename <- ""
        cytof_qc_file_statuses$unsuccessful_abnormal_gating_unflag_filename <- ""
        cytof_qc_file_statuses$successful_abnormal_gating_unflag_filename <- ""
        cytof_qc_file_statuses$unsuccessful_updated_qc_report_filename <- ""
        cytof_qc_file_statuses$successful_report_export_filenames <- ""
        cytof_qc_file_statuses$unsuccessful_report_export_filenames <- ""

        # We store a list of pre-processed fcs data and corresponding file names
        # for qc reports that were successfully exported in order to 
        # allow the user to flag gating as abnormal.
        cytof_qc_gating_inspection$pre_processed_data <- vector("list", 1000)
        cytof_qc_gating_inspection$currently_rendered_gating_filename <- ""
      }

      for (i in 1:num_files_uploaded) {
        fcs_data_path <- fcs_file_data_frame[i,]$datapath

        # We check to see if uploaded files have FCS file type. We do not 
        # check if "BCKG" is in filename since that is not always the case.
        if (!flowCore::isFCSfile(fcs_data_path)) {
          cytof_qc_control_var$is_uploaded_file_type_valid <- FALSE
          return()
        }
      }

      withProgress(value = 0, {
        for (i in 1:num_files_uploaded) {
          incProgress(amount = (1 / num_files_uploaded), message = paste("Generating",
          "QC report and exporting data to target directory for file", i, "/", 
          num_files_uploaded))
          fcs_filename <- fcs_file_data_frame[i,]$name
          fcs_data <- fcs_import_file_error_handler(fcs_file_data_frame[i,]$datapath)

          if (is.null(fcs_data)) {
            cytof_qc_control_var$fcs_file_import_error <- TRUE
            cytof_qc_file_statuses$unsuccessful_fcs_file_import_filenames <- c(cytof_qc_file_statuses$unsuccessful_fcs_file_import_filenames,
                                                                        fcs_filename)
            next
          }

          fcs_data_pre_processing <- mass_cytometry_pre_processing_error_handler(fcs_data$data, cofactor, 
                              event_channel, dna_channel, 
                              bead_gates)

          if (is.null(fcs_data_pre_processing)) {
            cytof_qc_control_var$pre_processing_error <- TRUE
            cytof_qc_file_statuses$unsuccessful_pre_processing_filenames <- c(cytof_qc_file_statuses$unsuccessful_pre_processing_filenames,
                                                                        fcs_filename)
            next
          }

          # This will be updated every time our pre-processing algorithm
          # is updated.
          qc_reporter_version <- "QCToolkit_v170622"

          qc_report <- generate_qc_report_error_handler(fcs_data, 
                                                        fcs_filename, 
                                                        fcs_data_pre_processing,
                                                        cofactor,
                                                        qc_reporter_version)

          if (is.null(qc_report)) {
            cytof_qc_control_var$qc_report_generation_error <- TRUE
            cytof_qc_file_statuses$unsuccessful_report_generation_filenames <- c(cytof_qc_file_statuses$unsuccessful_report_generation_filenames,
                                                                        fcs_filename)
            next
          }

          qc_report_export_status <- qc_report_export_error_handler(
                                                            cytof_qc_file_statuses$cytof_qc_report_dir,
                                                            qc_report,
                                                            fcs_filename)

          if (is.null(qc_report_export_status)) {
            cytof_qc_file_statuses$successful_report_export_filenames <- c(cytof_qc_file_statuses$successful_report_export_filenames,
                                                                                    fcs_filename)
            prepare_for_gating_inspection(fcs_filename, fcs_data_pre_processing, cytof_qc_gating_inspection)
            cytof_qc_gating_inspection$cytof_qc_report_tables[[fcs_filename]] <- qc_report
          } else if (sample_background_report_export_status == "Unsuccessful"){
            cytof_qc_control_var$qc_report_export_error <- TRUE
            cytof_qc_file_statuses$unsuccessful_report_export_filenames <- c(cytof_qc_file_statuses$unsuccessful_report_export_filenames,
                                                                                      fcs_filename)
          }    

          rm(fcs_data, 
            fcs_data_pre_processing, 
            qc_report)
        }
      })

      if (length(cytof_qc_file_statuses$successful_report_export_filenames) > 1) {
        cytof_qc_control_var$qc_report_export_success <- TRUE
        # We update our control variables to render the QC report table and 
        # gating inspection/visualization information on the UI for successfully 
        # completed files only.
        cytof_qc_control_var$render_gating_inspection <- TRUE
      }
    }
  })
}