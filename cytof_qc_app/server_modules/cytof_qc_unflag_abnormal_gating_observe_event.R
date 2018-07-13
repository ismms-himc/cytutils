unflag_abnormal_gating_observe_event <- function(input, cytof_qc_control_var, cytof_qc_file_statuses) {
	observeEvent(input$unflag_abnormal_gating, {
		cat('IN unflag_abnormal_gating_observe_event')
		# The name of the file we accidentally flagged as having abnormal gating
	    normal_gating_filename <- cytof_qc_file_statuses$successful_abnormal_gating_flag_filename
	    cytof_qc_report_dir <- cytof_qc_file_statuses$cytof_qc_report_dir

	    withProgress(
	    	message = "Undoing abnormal gating flag in previously exported QC report", {
    		      for (i in 1:7) {
			        incProgress(1/7)
			        Sys.sleep(0.25)
			      }
			    abnormal_gating_unflag_transfer_status <- unflag_abnormal_gating_in_previously_exported_qc_report_error_handler(normal_gating_filename,
			    																												cytof_qc_report_dir)
	    	}
    	)

    	cat(paste("abnormal_gating_unflag_transfer_status: ", abnormal_gating_unflag_transfer_status))

        if (abnormal_gating_unflag_transfer_status == "Unsuccessful") {
	        cytof_qc_control_var$abnormal_gating_unflag_error <- TRUE
        	cytof_qc_file_statuses$unsuccessful_abnormal_gating_unflag_filename <- normal_gating_filename
			# We clear old success messages.
        	cytof_qc_control_var$successful_abnormal_gating_unflag <- FALSE
			cytof_qc_file_statuses$successful_abnormal_gating_unflag_filename <- ""
        } else {
            cytof_qc_control_var$successful_abnormal_gating_unflag <- TRUE
        	cytof_qc_file_statuses$successful_abnormal_gating_unflag_filename <- normal_gating_filename
			# We clear old success messages.
        	cytof_qc_control_var$successful_abnormal_gating_flag <- FALSE
			cytof_qc_file_statuses$successful_abnormal_gating_flag_filename <- ""
			# We clear old error messages.
	        cytof_qc_control_var$abnormal_gating_unflag_error <- FALSE
        	cytof_qc_file_statuses$unsuccessful_abnormal_gating_unflag_filename <- ""
        }
	})
}
