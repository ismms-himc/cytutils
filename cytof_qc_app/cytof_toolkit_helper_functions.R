fcs_import_file <- function(filename,
                            return_original = FALSE,
                            transformation = TRUE) {

  fcs_data <- flowCore::read.FCS(filename, transformation)
  data <- dplyr::as_data_frame(as.data.frame(fcs_data@exprs))
  colnames(data) <- paste0(fcs_data@parameters@data$name, "_",
                           fcs_data@parameters@data$desc)
  # rename columns to be R-friendly
  colnames(data) <- gsub("-| |#", "_", colnames(data))
  
  if (return_original) {
    list(
      data = data,
      obj = fcs_data
    )
  } else {
    dplyr::as_data_frame(data)
  }
}

fcs_import_file_error_handler <- function(filename) {
  tryCatch(fcs_import_file(filename, 
                            return_original = TRUE),
    error = function(e) { NULL }
  )
}

mass_cytometry_pre_processing <- function(
  fcs_data, cofactor, event_channel, dna_channel, bead_gates,
  debris_num_gmms = 10, debris_gmm_subsample_n = 1000, debris_vote_thresh = 0.5,
  verbose = TRUE) {

    
  start_t <- proc.time()
 
  # find the explicit channel name for each mass
  fcs_event_channel <- fcs_find_channels(fcs_data, event_channel)
  fcs_dna_channel <- fcs_find_channels(fcs_data, dna_channel)
  fcs_bead_channels <- fcs_find_channels(fcs_data, bead_gates$channel)
  # transform the data and add an index for later tracking
  gating_data <- fcs_data[, c(fcs_event_channel,
                              fcs_dna_channel,
                              fcs_bead_channels)]
  gating_data[, c(fcs_dna_channel, fcs_bead_channels)] <-
    asinh(gating_data[, c(fcs_dna_channel, fcs_bead_channels)] / cofactor)
  gating_data$index <- 1:nrow(gating_data)
  
  # step 1: identify beads and bead-cell doublets using threshold
  if (verbose) cat("identifying beads\n")
  bead_data <- gating_data
  for (bead_gate_idx in 1:nrow(bead_gates)) {
    bead_gate <- bead_gates[bead_gate_idx, ]
    bead_gate_channel <- fcs_find_channels(gating_data, bead_gate$channel)
    bead_gate_data <- bead_data[[bead_gate_channel]]
    bead_gate_indices <- which(bead_gate$min < bead_gate_data &
                                 bead_gate_data < bead_gate$max)
    bead_data <- bead_data[bead_gate_indices, ]

  }
  
  bead_indices <- bead_data$index

  gating_data <- dplyr::filter(gating_data, !(index %in% bead_indices))
  # step 2: separate cells from debris
  if (verbose) cat("identifying debris ... ")
  t <- proc.time()
  cell_indices <- gating_data$index
  dna <- gating_data[[fcs_dna_channel]]
  min_dna_index <- head(which(dna == min(dna)), 1)
  
  debris_votes <- lapply(1:debris_num_gmms, function(iter) {
    # fit a GMM with three components to a subset of the data and use it to
    # assign all of the data. lowest cluster is classified as debris
    indices <- sample(1:length(dna), min(debris_gmm_subsample_n, length(dna)))
    gmm <- Mclust(dna[indices], G = 3)
    assignment <- predict(gmm, dna)$classification
    min_cluster <- assignment[min_dna_index]
    # sometimes the GMM finds a very "wide" cluster that covers debris and some
    # doublets. in order to avoid removing doublets, cells with DNA intensity
    # higher than most populated cluster are not removed
    assignment_counts <- table(assignment)
    most_populated_mean <- gmm$parameters$mean[
      which(assignment_counts == max(assignment_counts))]
    # return debris identification for each cell
    dplyr::data_frame(
      iter = iter,
      index = cell_indices,
      debris = (assignment == min_cluster) & (dna < most_populated_mean)
    )
  }) %>% dplyr::bind_rows()
  
  debris_indices <- debris_votes %>%
    dplyr::group_by(index) %>%
    dplyr::summarize(per_debris = mean(debris)) %>%
    dplyr::filter(per_debris >= debris_vote_thresh) %>%
    `$`(index)
  
  debris_data <- dplyr::filter(gating_data, index %in% debris_indices)
  gating_data <- dplyr::filter(gating_data, !(index %in% debris_indices))
  
  cat(paste0(round((proc.time() - t)[3]), " seconds\n"))

  # step 3: separate bead + cell doublets
  bead_singlet_indices <- which(bead_data[[fcs_dna_channel]] <
                                  min(gating_data[[fcs_dna_channel]]))

  bead_cell_doublet_indices <- which(bead_data[[fcs_dna_channel]] >
                                min(gating_data[[fcs_dna_channel]]))
  
  bead_cell_doublets_data <- bead_data[bead_cell_doublet_indices, ]
  
  bead_data <- bead_data[bead_singlet_indices, ]

  list(
    bead_data = bead_data,
    debris_data = debris_data,
    gating_data = gating_data,
    bead_cell_doublets_data = bead_cell_doublets_data,
    cell_indices = gating_data$index
  )
}

mass_cytometry_pre_processing_error_handler <- function(fcs_data, cofactor, 
                            event_channel, dna_channel, 
                            bead_gates) {
  tryCatch(mass_cytometry_pre_processing(fcs_data, cofactor, 
                            event_channel, dna_channel, 
                            bead_gates),
    error = function(e) { NULL }
  )
}

generate_sample_background_report <- function(background_fcs_data, 
                              background_fcs_data_filename, 
                              background_fcs_data_pre_processing, 
                              cofactor) {
  sample_background_report_timestamp <- Sys.time()

  gating_data <- background_fcs_data$data[background_fcs_data_pre_processing$gating_data$index, ]
  gating_data_channel_names <- colnames(gating_data)
  gating_data_channel_names <- tail(gating_data_channel_names, length(gating_data_channel_names) - 2)
  gating_data_channel_names <- gsub(".{5}$", "", gating_data_channel_names)
  
  # We manually fix the channel names for Ir191Di_DNA and Ir193Di_DNA to avoid
  # a trailing "D" in the channel name (i.e. Ir191D)
  gating_data_channel_names[58] <- "Ir191"
  gating_data_channel_names[60] <- "Ir193"

  num_gating_channels <- length(gating_data_channel_names)

  num_gating_channels_data_matrix_rows <- length(gating_data[[fcs_find_channels(gating_data, 
                                            "Y89")]])

  # Note: Pre-determining # of rows as per below assumes that each channel has the 
  # same number of entries.
  gating_channels_data_matrix <- matrix(nrow = num_gating_channels_data_matrix_rows, ncol = 66)
  colnames(gating_channels_data_matrix) <- gating_data_channel_names

  # Vector of data for each channel.
  for (i in 1:num_gating_channels) {
    current_channel_gating_data <- gating_data[[fcs_find_channels(gating_data, gating_data_channel_names[i])]]
    gating_channels_data_matrix[, i] <- current_channel_gating_data
  }

  gating_data_channel_medians <- matrixStats::colMedians(gating_channels_data_matrix)
  sample_background_report <-
    dplyr::data_frame(
      `Sample Background Report Timestamp` = sample_background_report_timestamp,
      Filename = background_fcs_data_filename,
      `Start Time` =
        paste0(background_fcs_data$obj@description$`$DATE`, " ",
               background_fcs_data$obj@description$`$BTIM`),
      `End Time` =
        paste0(background_fcs_data$obj@description$`$DATE`, " ",
               background_fcs_data$obj@description$`$ETIM`),
      `Total Events` = nrow(background_fcs_data$data),
      `Total Cells` = nrow(gating_data),
      `Y89 Median` = median(gating_channels_data_matrix[, 'Y89']),
      `Nb93 Median` = median(gating_channels_data_matrix[, 'Nb93']),
      `Pd102 Median` = median(gating_channels_data_matrix[, 'Pd102']),
      `Rh103 Median` = median(gating_channels_data_matrix[, 'Rh103']),
      `Pd104 Median` = median(gating_channels_data_matrix[, 'Pd104']),
      `Pd105 Median` = median(gating_channels_data_matrix[, 'Pd105']),
      `Pd106 Median` = median(gating_channels_data_matrix[, 'Pd106']),
      `Pd108 Median` = median(gating_channels_data_matrix[, 'Pd108']),
      `Pd110 Median` = median(gating_channels_data_matrix[, 'Pd110']),
      `In113 Median` = median(gating_channels_data_matrix[, 'In113']),
      `In115 Median` = median(gating_channels_data_matrix[, 'In115']),
      `Sn120 Median` = median(gating_channels_data_matrix[, 'Sn120']),
      `I127 Median` = median(gating_channels_data_matrix[, 'I127']),
      `Xe131 Median` = median(gating_channels_data_matrix[, 'Xe131']),
      `Cs133 Median` = median(gating_channels_data_matrix[, 'Cs133']),
      `Ba138 Median` = median(gating_channels_data_matrix[, 'Ba138']),
      `La139 Median` = median(gating_channels_data_matrix[, 'La139']),
      `Ce140 Median` = median(gating_channels_data_matrix[, 'Ce140']),
      `Pr141 Median` = median(gating_channels_data_matrix[, 'Pr141']),
      `Nd142 Median` = median(gating_channels_data_matrix[, 'Nd142']),
      `Ce142 Median` = median(gating_channels_data_matrix[, 'Ce142']),
      `Nd143 Median` = median(gating_channels_data_matrix[, 'Nd143']),
      `Nd144 Median` = median(gating_channels_data_matrix[, 'Nd144']),
      `Nd145 Median` = median(gating_channels_data_matrix[, 'Nd145']),
      `Nd146 Median` = median(gating_channels_data_matrix[, 'Nd146']),
      `Sm147 Median` = median(gating_channels_data_matrix[, 'Sm147']),
      `Nd148 Median` = median(gating_channels_data_matrix[, 'Nd148']),
      `Sm149 Median` = median(gating_channels_data_matrix[, 'Sm149']),
      `Nd150 Median` = median(gating_channels_data_matrix[, 'Nd150']),
      `Eu151 Median` = median(gating_channels_data_matrix[, 'Eu151']),
      `Sm152 Median` = median(gating_channels_data_matrix[, 'Sm152']),
      `Eu153 Median` = median(gating_channels_data_matrix[, 'Eu153']),
      `Sm154 Median` = median(gating_channels_data_matrix[, 'Sm154']),
      `Gd155 Median` = median(gating_channels_data_matrix[, 'Gd155']),
      `Gd156 Median` = median(gating_channels_data_matrix[, 'Gd156']),
      `Gd158 Median` = median(gating_channels_data_matrix[, 'Gd158']),
      `Tb159 Median` = median(gating_channels_data_matrix[, 'Tb159']),
      `Gd160 Median` = median(gating_channels_data_matrix[, 'Gd160']),
      `Dy161 Median` = median(gating_channels_data_matrix[, 'Dy161']),
      `Dy162 Median` = median(gating_channels_data_matrix[, 'Dy162']),
      `Dy163 Median` = median(gating_channels_data_matrix[, 'Dy163']),
      `Dy164 Median` = median(gating_channels_data_matrix[, 'Dy164']),
      `Ho165 Median` = median(gating_channels_data_matrix[, 'Ho165']),
      `Er166 Median` = median(gating_channels_data_matrix[, 'Er166']),
      `Er167 Median` = median(gating_channels_data_matrix[, 'Er167']),
      `Er168 Median` = median(gating_channels_data_matrix[, 'Er168']),
      `Tm169 Median` = median(gating_channels_data_matrix[, 'Tm169']),
      `Er170 Median` = median(gating_channels_data_matrix[, 'Er170']),
      `Yb171 Median` = median(gating_channels_data_matrix[, 'Yb171']),
      `Yb172 Median` = median(gating_channels_data_matrix[, 'Yb172']),
      `Yb173 Median` = median(gating_channels_data_matrix[, 'Yb173']),
      `Yb174 Median` = median(gating_channels_data_matrix[, 'Yb174']),
      `Lu175 Median` = median(gating_channels_data_matrix[, 'Lu175']),
      `Lu176 Median` = median(gating_channels_data_matrix[, 'Lu176']),
      `Yb176 Median` = median(gating_channels_data_matrix[, 'Yb176']),
      `Ta181 Median` = median(gating_channels_data_matrix[, 'Ta181']),
      `Os189 Median` = median(gating_channels_data_matrix[, 'Os189']),
      `Ir191 Median` = median(gating_channels_data_matrix[, 'Ir191']),
      `Os192 Median` = median(gating_channels_data_matrix[, 'Os192']),
      `Ir193 Median` = median(gating_channels_data_matrix[, 'Ir193']),
      `Pt194 Median` = median(gating_channels_data_matrix[, 'Pt194']),
      `Pt195 Median` = median(gating_channels_data_matrix[, 'Pt195']),
      `Pt196 Median` = median(gating_channels_data_matrix[, 'Pt196']),
      `Pt198 Median` = median(gating_channels_data_matrix[, 'Pt198']),
      `Pb208 Median` = median(gating_channels_data_matrix[, 'Pb208']),
      `Bi209 Median` = median(gating_channels_data_matrix[, 'Bi209']),
      `Sum of Medians` = sum(gating_data_channel_medians)
    )
}

generate_sample_background_report_error_handler <- function(background_fcs_data, 
                              background_fcs_data_filename, 
                              background_fcs_data_pre_processing, 
                              cofactor) {
  tryCatch(generate_sample_background_report(background_fcs_data, 
                        background_fcs_data_filename, 
                        background_fcs_data_pre_processing,
                        cofactor),
    error = function(e) { NULL }
  )
}


sample_background_report_export_error_handler <- function(sample_background_report_dir_path,
                                                            sample_background_report_all) {
  filepath <- paste0(sample_background_report_dir_path, "/sample_background_report.csv")
  tryCatch(write.csv(sample_background_report_all, 
                      filepath),
    error = function(e) { "Unsuccessful" })
}


# http://jeromyanglim.tumblr.com/post/33418725712/how-to-source-all-r-files-in-a-directory
source_dir <- function(path, trace = TRUE, ...) {
    for (nm in list.files(path, pattern = "\\.[RrSsQq]$")) {
       if(trace) cat(nm,":")           
       source(file.path(path, nm), ...)
       if(trace) cat("\n")
    }
}
