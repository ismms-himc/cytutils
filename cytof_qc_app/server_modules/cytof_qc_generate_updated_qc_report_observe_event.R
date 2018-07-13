# default gating configuration
event_channel <- "Event_length"
dna_channel <- "191Di"
bead_gates <- dplyr::data_frame(
  channel = c("140Di", "175Di"),
  min = c(5, 5),
  max = c(Inf, Inf)
)
cofactor <- 5

cytof_qc_generate_updated_qc_report_observe_event <- function(input, output, cytof_qc_control_var, cytof_qc_file_statuses, cytof_qc_gating_inspection) {

	observeEvent(input$generate_updated_qc_report, {

		# We remove old updated QC report generation/export success
		# and error messages
	    cytof_qc_control_var$updated_qc_report_generation_error <- FALSE
		cytof_qc_file_statuses$unsuccessful_updated_qc_report_filename <- ""
	    cytof_qc_control_var$updated_qc_report_export_error <- FALSE
	    cytof_qc_file_statuses$unsuccessful_updated_qc_report_export_filename <- ""
	    cytof_qc_control_var$updated_qc_report_export_success <- FALSE
	    cytof_qc_file_statuses$successful_updated_qc_report_export_filename <- ""

		fcs_filename <- cytof_qc_gating_inspection$currently_rendered_gating_filename

	    fcs_file_data_frame <- input$cytof_qc_file
	    fcs_file_data_row <- which(fcs_file_data_frame$name == fcs_filename)
	    # We do not keep track of fcs_import_file errors or filenames associated
	    # with successful and unsuccessful fcs file imports at this point because
	    # the application will not render a gating visualization for bad files.
	    # Thus, we can safely assume there will not be a fcs file import error.
	    fcs_data <- fcs_import_file_error_handler(fcs_file_data_frame[fcs_file_data_row,]$datapath)

	    searchable_data_frames <- cytof_qc_gating_inspection$pre_processed_data

	    for (i in 1:length(searchable_data_frames)) {
	    	if (searchable_data_frames[[i]]$filename == fcs_filename) {
			    pre_processed_manually_gated_data <- searchable_data_frames[[i]]
			    break
	    	}
	    }

	    # We format our pre-processed, manually gated data so that it looks more
	    # similar to the output of the himc::mass_cytometry_pre_processing function.
	    fcs_data_pre_processing <- reformat_manually_gated_pre_processed_data(pre_processed_manually_gated_data)

        qc_reporter_version <- "QCToolkit_v170622"

        updated_qc_report <- generate_qc_report_error_handler(fcs_data, 
                                                      fcs_filename, 
                                                      fcs_data_pre_processing,
                                                      cofactor,
                                                      qc_reporter_version,
                                                      ab_gate = "True")

        if (is.null(updated_qc_report)) {
          cytof_qc_control_var$updated_qc_report_generation_error <- TRUE
          cytof_qc_file_statuses$unsuccessful_updated_qc_report_filename <- fcs_filename
          return()
        }

	    withProgress(
	    	message = "Generating and exporting updated QC report", {
	    		      for (i in 1:7) {
				        incProgress(1/7)
				        Sys.sleep(0.25)
				      }

	        qc_report_export_status <- qc_report_export_error_handler(
                                                            cytof_qc_file_statuses$cytof_qc_report_dir,
                                                            updated_qc_report,
                                                            fcs_filename)
	    	}
    	)

	    if (is.null(qc_report_export_status)) {
          cytof_qc_control_var$updated_qc_report_export_success <- TRUE
          cytof_qc_file_statuses$successful_updated_qc_report_export_filename <- fcs_filename
	    } else if (sample_background_report_export_status == "Unsuccessful"){
          cytof_qc_control_var$updated_qc_report_export_error <- TRUE
          cytof_qc_file_statuses$unsuccessful_updated_qc_report_export_filename <- fcs_filename
        }    

        rm(fcs_data, 
          fcs_data_pre_processing, 
          updated_qc_report,
          qc_report_export_status)
	})
}