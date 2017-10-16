# channel_population_relationships_filepath <- "channel_population_relationships.csv"
# base_fcs_data_filepath <- "170203_IOF_1_1_singlets.fcs"

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
		neg_indices <- set_diff(neg_indices, pos_indices)

		aof_for_current_channel <- calculateAof(x, pos_indices, neg_indices)
		# print(paste("aof_for_current_channel:", aof_for_current_channel))
		aof_results_row <- data.frame(channel, aof_for_current_channel)
		names(aof_results_row) <- c("Channel Name", "Aof")
		# print(paste("aof_results_row:", aof_results_row))

		aof_results <- rbind(aof_results, aof_results_row)
		# View(aof_results)

		i <- i + 1
	}

	aof_results
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