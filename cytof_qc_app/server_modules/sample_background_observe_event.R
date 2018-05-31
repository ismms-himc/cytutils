library(googlesheets)
library(flowCore)
library(himc)
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
shinyDirChoose(input, id = "sample_background_report_dir", roots = c(home = '~'))

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
      sample_background_control_var$aggr_sample_background_report_exists <- TRUE
      sample_background_control_var$aggr_sample_background_report_export_error <- FALSE
      sample_background_control_var$aggr_sample_background_report_export_success <- FALSE


      # We also reset the reactive values of our sample_background_file_statuses variable to
      # ensure error messages include updated file names after a user
      # re-uploads files.
      sample_background_file_statuses$unsuccessful_pre_processing_filenames <- ""
      sample_background_file_statuses$unsuccessful_report_generation_filenames <- ""
      sample_background_file_statuses$unsuccessful_fcs_file_import_filenames <- ""
      sample_background_file_statuses$successful_sample_background_report_completion_filenames <- ""
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

      sample_background_report_all <- data.frame(matrix(ncol = 73, nrow = 0))
      cols <- c("Sample Background Report Timestamp", 
                "Filename", 
                "Start Time",
                "End Time",
                "Total Events",
                "Total Cells",
                "Y89 Median",
                "Nb93 Median",
                "Pd102 Median",
                "Rh103 Median",
                "Pd104 Median",
                "Pd105 Median",
                "Pd106 Median",
                "Pd108 Median",
                "Pd110 Median",
                "In113 Median",
                "In115 Median",
                "Sn120 Median",
                "I127 Median",
                "Xe131 Median",
                "Cs133 Median",
                "Ba138 Median",
                "Ce140 Median",
                "La139 Median",
                "Pr141 Median",
                "Nd142 Median",
                "Ce142 Median",
                "Nd143 Median",
                "Nd144 Median",
                "Nd145 Median",
                "Nd146 Median",
                "Sm147 Median",
                "Nd148 Median",
                "Sm149 Median",
                "Nd150 Median",
                "Eu151 Median",
                "Sm152 Median",
                "Eu153 Median",
                "Sm154 Median",
                "Gd155 Median",
                "Gd156 Median",
                "Gd158 Median",
                "Tb159 Median",
                "Gd160 Median",
                "Dy161 Median",
                "Dy162 Median",
                "Dy163 Median",
                "Dy164 Median",
                "Ho165 Median",
                "Er166 Median",
                "Er167 Median",
                "Er168 Median",
                "Tm169 Median",
                "Er170 Median",
                "Yb171 Median",
                "Yb172 Median",
                "Yb173 Median",
                "Yb174 Median",
                "Lu175 Median",
                "Lu176 Median",
                "Yb176 Median",
                "Ta181 Median",
                "Os189 Median",
                "Ir191 Median",
                "Os192 Median",
                "Ir193 Median",
                "Pt194 Median",
                "Pt195 Median",
                "Pt196 Median",
                "Pt198 Median",
                "Pb208 Median",
                "Bi209 Median",
                "Sum of Medians"
              )
      colnames(sample_background_report_all) <- cols

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
        }

        sample_background_report_all <- rbind(sample_background_report_all, sample_background_report)
        sample_background_file_statuses$successful_sample_background_report_completion_filenames <- c(sample_background_file_statuses$successful_sample_background_report_completion_filenames, 
                                                                  background_data_filename)

        rm(background_fcs_data_pre_processing,
            sample_background_report)
      }


      if (nrow(sample_background_report_all) > 0) {
        sample_background_report_export_status <- sample_background_report_export_error_handler(
                                                            sample_background_file_statuses$sample_background_report_dir,
                                                            sample_background_report_all)

        if (is.null(sample_background_report_export_status)) {
          sample_background_control_var$aggr_sample_background_report_export_success <- TRUE

        } else if (sample_background_report_export_status == "Unsuccessful"){
          sample_background_control_var$aggr_sample_background_report_export_error <- TRUE
        }
      } else {
        # A background report was not successfully generated for any of our samples
        # and thus, there is no report to export.
        sample_background_control_var$aggr_sample_background_report_exists <- FALSE
      }
    })
  })
}