channel_population_relationships_filepath <- "channel_population_relationships.csv"
base_fcs_data_filepath <- "170203_IOF_1_1_singlets.fcs"
sample_id <- "iof_1"
# Use manually gated data to assign sample cells to specific populations
sample_1_base_fcs_data <- flowCore::read.FCS("170203_IOF_1_1_singlets.fcs")
manual_labeling_filepath <- "LabelCellsByManualGatingTest.manual_labeling.csv"
samples_filepath <- "LabelCellsByManualGatingTest.csv"
data_dir <- "/home/himc/himc/cytutils-test-files/test6"

single_sample_labels <- generatePopulationAssignments(
							manual_labeling_filepath, 
							samples_filepath, 
							data_dir)

# TODO: fn name ok?
# TODO: test calculated aof when csv Er168Di TRUE is t_cell only
# NB: population names should be the same as in our manual labeling csv
calculateMultiChannelAof <- function(channel_population_relationships_filepath, base_fcs_data_filepath, single_sample_labels, sample_id) {

	channel_population_relationships <- .readChannelPopulationRelationships(channel_population_relationships_filepath)
	channels <- channel_population_relationships$channel
	num_populations <- length(colnames(channel_population_relationships))
	# We remove the first column, "channel" from population names.
	populations <- colnames(channel_population_relationships)[2:num_populations]
	base_fcs_data <- flowCore::read.FCS(base_fcs_data_filepath)
	y <- single_sample_labels[sample_id]
	y <- data.frame(y, stringsAsFactors = FALSE)

	aof_results <- data.frame(matrix(ncol = 2, nrow = 0), stringsAsFactors = FALSE)
	aof_results_colnames <- c("Channel Name", "Aof")
	colnames(aof_results) <- aof_results_colnames

	i <- 1
	for (channel in channels) {
		x <- base_fcs_data@exprs[, channel]
		# We subtract one below to account for the fact that the first column is
		# "channel" and not a population.
		# TODO: maybe change below var names to pos_population_indices instead of channel_population
		# print(paste("i", i))
		print(paste("channel:", channel))

		pos_channel_population_indices <- grep(TRUE, channel_population_relationships[i,]) - 1
		# print(paste("pos_channel_population_indices", pos_channel_population_indices))

		pos_channel_populations <- populations[pos_channel_population_indices]
		# print(paste("pos_channel_populations", pos_channel_populations))

		neg_channel_population_indices <- grep(FALSE, channel_population_relationships[i,]) - 1
		# print(paste("neg_channel_population_indices", neg_channel_population_indices))

		neg_channel_populations <- populations[neg_channel_population_indices]
		# print(paste("neg_channel_populations", neg_channel_populations))

		
		pos_indices <- c()
		# print(paste("pos_indices", pos_indices))

		for (pos_channel_population in pos_channel_populations) {
			# TODO: check for scoping issues
			# print(paste("pos_channel_population", pos_channel_population))
			target_col_name <- paste(sample_id, pos_channel_population, sep = ".")
			target_col_idx <- grep(target_col_name, colnames(y))
			# print(paste("y$sample_id$pos_channel_population", y$sample_id$pos_channel_population))
			pos_indices <- c(pos_indices, grep(TRUE, y[target_col_idx][,1]))
			# print(paste("pos_indices inside for", pos_indices))
		}

		pos_indices <- unique(pos_indices)
		# print(paste("pos_indices after for", pos_indices))

		neg_indices <- c()
		# print(paste("neg_indices", neg_indices))

		for (neg_channel_population in neg_channel_populations) {
			# print(paste("neg_channel_population", neg_channel_population))
			# print(paste("y$sample_id$neg_channel_population", y$sample_id$neg_channel_population))
			target_col_name <- paste(sample_id, neg_channel_population, sep = ".")
			target_col_idx <- grep(target_col_name, colnames(y))

			# TODO: check for scoping issues
			neg_indices <- c(neg_indices, grep(TRUE, y[target_col_idx][,1]))
			# print(paste("neg_indices inside for", neg_indices))

		}

		neg_indices <- unique(neg_indices)
		# We remove indices that are in our pos_indices vector.
		neg_indices <- setdiff(neg_indices, pos_indices)
		# print(paste("neg_indices after for", neg_indices))
		

		aof_for_current_channel <- calculateAof(x, pos_indices, neg_indices)
		print(paste("aof_for_current_channel:", aof_for_current_channel))
		aof_results_row <- data.frame(channel, aof_for_current_channel)
		names(aof_results_row) <- c("Channel Name", "Aof")
		# print(paste("aof_results_row:", aof_results_row))

		aof_results <- rbind(aof_results, aof_results_row)
		# View(aof_results)

		i <- i + 1
	}

	aof_results
}



calculateAof <- function(x, pos_indices, neg_indices, width = 0.05, cofactor = NULL) {
  if (length(pos_indices) == 0) {
    stop("no cells in positive population")
  }
  if (length(pos_indices) == length(x)) {
    stop("no cells in negative population")
  }
  if (max(pos_indices) > length(x)) {
    stop("positive population indices include values higher than length of x")
  }

  if (!is.null(cofactor)) x <- asinh(x / cofactor)

  # Set up positive and negative populations.
  pos <- x[pos_indices]
  neg <- x[neg_indices]

  if (mean(pos) < mean(neg)) {
    stop("positive population mean is lower than negative population mean")
  }

  # Calculate AOF thresholds.
  if (length(neg) == 1) {
    neg_high <- neg
  } else {
    neg_high <- qnorm(1 - width, mean(neg), sd(neg))
  }

  if (length(pos) == 1) {
    pos_low <- pos
  } else {
    pos_low <- qnorm(width, mean(pos), sd(pos))
  }

  # Calculate AOF.
  mean(c(mean(pos <= neg_high), mean(neg >= pos_low)))
}

.fcs_import_file <- function(filepath,
                            return_original = FALSE,
                            transformation = TRUE) {
  fcs_data <- flowCore::read.FCS(filepath, transformation)
  data <- dplyr::as_data_frame(as.data.frame(fcs_data@exprs))
  colnames(data) <- paste0(fcs_data@parameters@data$name, "_",
                           fcs_data@parameters@data$desc)
  # Rename columns to be R-friendly
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



generatePopulationAssignments <- function(manual_labeling_filepath, samples_filepath, data_dir) {
  manual_labeling <- read.csv(manual_labeling_filepath, stringsAsFactors = FALSE)
  labels <- setdiff(names(manual_labeling), c("sample_id", "base"))
  single_sample_analysis_labels <- list()
  samples <- .read_samples(samples_filepath)
  sample_ids <- samples$sample_id

  # Match manual gating labels to each sample.
  single_sample_labels <- lapply(sample_ids, function(sample_id) {
    cat(paste0(sample_id, "\n"))
    
    sample_labeling <- manual_labeling[manual_labeling$sample_id == sample_id, ]
    sample_base <- sample_labeling$base
    
    # Generate vector of clustering channels from base FCS file.
    base_fcs_data <- .fcs_import_file(file.path(data_dir, sample_base))
    base_fcs_data_cols <- colnames(base_fcs_data)
    time_col_index <- grep("Time", base_fcs_data_cols)
    time_col_name <- base_fcs_data_cols[time_col_index]
    event_length_col_index <- grep("Event_length", base_fcs_data_cols)
    event_length_col_name <- base_fcs_data_cols[event_length_col_index]
    non_channel_cols <- c(time_col_name, event_length_col_name)
    clustering_channels <- setdiff(base_fcs_data_cols, non_channel_cols)

    ssl <- lapply(labels, function(label) {
      cat(paste0("\t", label, "\n"))
      # Import FCS file for this label.
      label_fcs_data <-
        .fcs_import_file(file.path(data_dir, sample_labeling[[label]]))
      label_fcs_data[[label]] <- TRUE
      label_fcs_data_match <- 
        dplyr::left_join(base_fcs_data, label_fcs_data, by = clustering_channels)
      label_fcs_data_match[, label]
    }) %>% dplyr::bind_cols()
    
    ssl[is.na(ssl)] <- FALSE
    
    # ssl is a data frame with column names representing different labels (cell 
    # populations)
    ssl
  })

  names(single_sample_labels) <- sample_ids

  single_sample_labels
}


.greedyChannelAof <- function(x, y_indices, width = 0.05, cofactor = NULL) {
  # Cumulative frequency thresholds for search.
  search_thresh <- c(
    seq(0.2, 0.5, 0.05),
    seq(0.5, 0.8, 0.03),
    seq(0.8, 0.98, 0.02), 0.99
  )

  # Calculate mean and frequency of each cluster.
  y_labels <- names(y_indices)
  cluster_stats <- lapply(y_labels, function(y_label) {
    y_x <- x[y_indices[[y_label]]]
    data.frame(
      YLabel = y_label,
      Mean = mean(y_x),
      Freq = length(y_x) / length(x)
    )
  })
  cluster_stats <- do.call(rbind, cluster_stats)

  # Order clusters by mean x and calculate cumulative frequency.
  cluster_stats <- cluster_stats[order(cluster_stats$Mean), ]
  cluster_stats$CumFreq <- cumsum(cluster_stats$Freq)

  # Calculate AOF for each of the cumulative frequency thresholds.
  aof_values <- lapply(search_thresh, function(thresh) {
    # Find indices for all clusters positive for this threshold.
    pos_labels <-
      as.character(cluster_stats$YLabel[cluster_stats$CumFreq > thresh])
    if ((length(pos_labels) == 0) |
        (length(pos_labels) == nrow(cluster_stats))) {
      # Either positive or negative populations are empty.
      aof <- Inf
    } else {
      pos_indices <- unlist(y_indices[pos_labels])
      aof <- calculateAof(x, pos_indices, width = width, cofactor = cofactor)
    }

    aof
  })

  # Return the minimal AOF found.
  min(unlist(aof_values))
}


greedyCytometryAof <- function(fcs_data,
                               y,
                               channel_names = colnames(fcs_data),
                               width = 0.05,
                               cofactor = NULL,
                               verbose = TRUE) {
  if (length(y) != nrow(fcs_data)) {
    stop("y length should be equal to number of rows in fcs_data")
  }

  if (!all(channel_names %in% colnames(fcs_data))) {
    stop("channel_names should all be columns of fcs_data")
  }

  if (anyNA(fcs_data[, channel_names])) {
    stop("cannot calculate AOF when data includes missing values")
  }

  if (!is.factor(y)) {
    y <- as.factor(y)
  }

  # Find indices of cells in each cluster.
  if (verbose) message("Converting cluster assignments to indices vectors")
  y_labels <- unique(y)
  print(paste("y_labels:", y_labels))
  y_indices <- lapply(y_labels, function(y_label) {
    which(y == y_label)
  })
  names(y_indices) <- y_labels

  # Run greedy algorithm over each channel.
  aof_values <- lapply(channel_names, function(channel_name) {
    if (verbose) message(channel_name)
    x <- fcs_data[, channel_name]

    data.frame(
      ChannelName = channel_name,
      Aof = .greedyChannelAof(x, y_indices, width = width, cofactor = cofactor)
    )
  })

  do.call(rbind, aof_values)
}


.readChannelPopulationRelationships <- function(channel_population_relationships_filepath) {

  if (!file.exists(channel_population_relationships_filepath)) {
    stop(paste0("unable to find file ", channel_population_relationships_filepath))
  }
  
  channel_population_relationships <- read.csv(channel_population_relationships_filepath,  stringsAsFactors = FALSE)
  
  if (colnames(channel_population_relationships)[1] != "channel") {
    stop("first column of csv should be channel")
  }
  
  if (any(duplicated(channel_population_relationships$channel))) {
    stop("channel should be unique")
  }
  
  if (ncol(channel_population_relationships) == 1) {
    stop("samples file includes no cell population columns")
  }
  
  channel_population_relationships$channel <- as.character(channel_population_relationships$channel)
  
  channel_population_relationships
}