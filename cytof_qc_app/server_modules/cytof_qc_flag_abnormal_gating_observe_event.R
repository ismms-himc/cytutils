flag_abnormal_gating_observe_event <- function(input, cytof_qc_gating_inspection, cytof_qc_control_var, cytof_qc_file_statuses) {
	observeEvent(input$flag_abnormal_gating, {

		# We remove old updated QC report generation/export success
		# and error messages
	    cytof_qc_control_var$updated_qc_report_generation_error <- FALSE
		cytof_qc_file_statuses$unsuccessful_updated_qc_report_filename <- ""
	    cytof_qc_control_var$updated_qc_report_export_error <- FALSE
	    cytof_qc_file_statuses$unsuccessful_updated_qc_report_export_filename <- ""
	    cytof_qc_control_var$updated_qc_report_export_success <- FALSE
	    cytof_qc_file_statuses$successful_updated_qc_report_export_filename <- ""

	    abnormal_gating_filename <- cytof_qc_gating_inspection$currently_rendered_gating_filename
	    cytof_qc_report_dir <- cytof_qc_file_statuses$cytof_qc_report_dir
	    withProgress(
	    	message = "Flagging abnormal gating in previously exported QC report", {
	    		      for (i in 1:7) {
				        incProgress(1/7)
				        Sys.sleep(0.25)
				      }
			    abnormal_gating_flag_transfer_status <- flag_abnormal_gating_in_exported_qc_report_error_handler(abnormal_gating_filename,
			    																								cytof_qc_report_dir)
	    	}
    	)

        if (abnormal_gating_flag_transfer_status == "Unsuccessful") {
	        cytof_qc_control_var$abnormal_gating_flag_error <- TRUE
        	cytof_qc_file_statuses$unsuccessful_abnormal_gating_flag_filename <- abnormal_gating_filename
			# We clear old success messages.
        	cytof_qc_control_var$successful_abnormal_gating_flag <- FALSE
			cytof_qc_file_statuses$successful_abnormal_gating_flag_filename <- ""
        	cytof_qc_control_var$successful_abnormal_gating_unflag <- FALSE
			cytof_qc_file_statuses$successful_abnormal_gating_unflag_filename <- ""
        } else {
            cytof_qc_control_var$successful_abnormal_gating_flag <- TRUE
        	cytof_qc_file_statuses$successful_abnormal_gating_flag_filename <- abnormal_gating_filename
        	# We clear old success messages.
        	cytof_qc_control_var$successful_abnormal_gating_unflag <- FALSE
			cytof_qc_file_statuses$successful_abnormal_gating_unflag_filename <- ""
			# We clear old error messages.
	        cytof_qc_control_var$abnormal_gating_flag_error <- FALSE
        	cytof_qc_file_statuses$unsuccessful_abnormal_gating_flag_filename <- ""
        	# We render a section to allow users to manually gate data if they
        	# successfully flagged gating as abnormal.
        	cytof_qc_control_var$render_manual_gating <- TRUE
        }
	})
}