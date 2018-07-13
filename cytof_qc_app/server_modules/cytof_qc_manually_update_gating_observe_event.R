cytof_qc_manually_update_gating_observe_event <- function(input, output, session, cytof_qc_control_var, cytof_qc_file_statuses, cytof_qc_gating_inspection) {

	observeEvent(input$manually_update_gating, {
		# We remove old updated QC report generation/Google Drive transfer success
		# and error messages
	    cytof_qc_control_var$updated_qc_report_generation_error <- FALSE
		cytof_qc_file_statuses$unsuccessful_updated_qc_report_filename <- ""
	    cytof_qc_control_var$google_drive_updated_qc_report_transfer_error <- FALSE
	    cytof_qc_file_statuses$unsuccessful_google_drive_updated_qc_report_transfer_filename <- ""
	    cytof_qc_control_var$google_drive_updated_qc_report_transfer_success <- FALSE
	    cytof_qc_file_statuses$successful_google_drive_updated_qc_report_transfer_filename <- ""

	    gating_filename <- cytof_qc_gating_inspection$currently_rendered_gating_filename

	    # We look for the data frame corresponding to the current gating
	    # visualization.
	    searchable_data_frames <- cytof_qc_gating_inspection$pre_processed_data

	    for (i in 1:length(searchable_data_frames)) {
	    	if (searchable_data_frames[[i]]$filename == gating_filename) {
			    gating_visualization_data_frame <- searchable_data_frames[[i]]
			    gating_visualization_data_frame_index <- i
			    break
	    	}
	    }

	    if (exists("gating_visualization_data_frame")) {
	    	gated_population <- brushedPoints(gating_visualization_data_frame, input$manual_gating_brush)
	    	target_category <- input$manual_gating_category
	    	target_indices <- gated_population$index
	    	target_rows <- which(gating_visualization_data_frame$index %in% target_indices)

	    	# This conditional prevents an error-triggered premature exit from 
	    	# the application.
	    	if (is.integer(target_rows) & length(target_rows) == 0) {
	    		cytof_qc_control_var$manual_gating_error <- TRUE
				cytof_qc_control_var$manual_gating_success <- FALSE
	    	} else {
		    	gating_visualization_data_frame[target_rows,]$category <- target_category
		    	# We update the data frame held in cytof_qc_gating_inspection$pre_processed_data
		    	cytof_qc_gating_inspection$pre_processed_data[[gating_visualization_data_frame_index]] <- gating_visualization_data_frame
	    		cytof_qc_control_var$manual_gating_error <- FALSE
				cytof_qc_control_var$manual_gating_success <- TRUE
	    	}
	    }

	    output$gating_visualization <- renderPlot({
	      
	      	# We create a Progress object instead of using withProgress because
	      	# using the latter resulted in a progress bar that disappeared
	      	# too quickly.
	      	progress <- Progress$new(session, min=1, max=5)
		    on.exit(progress$close())

		    progress$set(message = 'Updating gating visualization')

		    gating_plot <- ggplot(data = gating_visualization_data_frame,
		    mapping = aes(x = Ir191Di_DNA,
		            y = Ce140Di_NA)) + geom_point(aes(colour = category)) + 
		  					ggtitle(gating_filename) + 
							  theme(plot.title = element_text(hjust = 0.5)) +
		                      xlab("Ir191Di DNA") +
		                      ylab("Ce140Di")

	        for (i in 1:5) {
		      progress$set(value = i)
		      Sys.sleep(0.5)
		    }

		    gating_plot
		})
	})
}