cytof_qc_gating_visualization_observe_event <- function(input, output, session, cytof_qc_control_var, cytof_qc_gating_inspection) {

	observeEvent(input$generate_gating_visualization, {
	    gating_filename <- input$gating_filename
	    if (gating_filename == "") {
	    	# TODO: Handle error. User must select filename from dropdown prior
	    	# to generating visualization.
	    	return
	    }

	    # We reset this to false in order to avoid a user toggling between data
	    # files and mistakenly associating an update QC report with an incorrect
	    # filename and visualization. 
	    cytof_qc_control_var$render_manual_gating <- FALSE
	    # We reset the two variables below in order to clear both the old success/
	    # error messages as well as the "Undo Abnormal Gating Flag" button.
	    cytof_qc_control_var$successful_abnormal_gating_flag <- FALSE
	    cytof_qc_control_var$abnormal_gating_flag_error <- FALSE
	    cytof_qc_control_var$manual_gating_success <- FALSE
	    cytof_qc_control_var$manual_gating_error <- FALSE

	    # We look for the data frame corresponding to the file the user chose to
	    # visualize via the drop down menu.
	    searchable_data_frames <- cytof_qc_gating_inspection$pre_processed_data

	    for (i in 1:length(searchable_data_frames)) {
	    	cat("length")
	    	cat(length(searchable_data_frames))
	    	if (searchable_data_frames[[i]]$filename == gating_filename) {
			    target_gating_visualization_data_frame <- searchable_data_frames[[i]]
			    break
	    	}
	    }

	    View(target_gating_visualization_data_frame)

	    if (exists("target_gating_visualization_data_frame")) {
		    output$gating_visualization <- renderPlot({
	  		    withProgress(
			    	message = "Generating gating visualization", {

					target_plot <- ggplot(data = target_gating_visualization_data_frame,
							    		mapping = aes(x = Ir191Di_DNA,
							            y = Ce140Di_NA)) + geom_point(aes(colour = category)) + 
							  					ggtitle(gating_filename) + 
												  theme(plot.title = element_text(hjust = 0.5)) +
							                      xlab("Ir191Di DNA") +
							                      ylab("Ce140Di")

			      	for (i in 1:5) {
				        incProgress(1/5)
				        Sys.sleep(0.25)
				   	}

					target_plot
		    	})
			})



			output$abnormal_gating_flag <- renderUI({
				div(
			        tags$hr(style = "border-color: #C0C0C0; margin: 35px 0px 27px 0px"),
					h4("Abnormal Gating?"),
			        actionButton(inputId = "flag_abnormal_gating",
		                label = "Flag Abnormal Gating", 
		                icon = icon("flag"),
		                style = "margin-top: 7px;")
					)
			})

			cytof_qc_gating_inspection$currently_rendered_gating_filename <- gating_filename
	    }
	})
}