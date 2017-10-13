#' Average Overlap Frequency (AOF).
#'
#' Calculate the Average Overlap Frequency (AOF) for a single cytometry channel.
#'
#' @param x Numeric vector corresponding to one channel in cytometry data.
#' @param pos_indices Indices of cells positive for this channel.
#' @param neg_indices Indices of cells negative for this channel.
#' @param width Width of high threshold of negative population and low
#' threshold of positive population.
#' @param cofactor If supplied, data will be transformed using inverse
#' hyperbolic sin with given cofactor.
#' @return The AOF between positive and negative populations for x.
#' @export
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


#' Read samples file
#'
#' Read the samples CSV file. The file includes a row for every sample in the
#' experiment. The first column is always sample_id, a unique sample ID. The
#' second column, filename, contains the name of the "base" FCS file (i.e. the
#' FCS file with data for all cells from the sample, regardless of manual gating
#' population designation).

#' Calculate the Average Overlap Frequency (AOF) for a single cytometry channel.
#'
#' @param samples_filepath Filepath of the samples csv file that outlines sample 
#' IDs and corresponding FCS file paths.
#' @return A data frame with the first column representing the sample_id and the
#' second column representing the base FCS filename corresponding to that sample_id.
.read_samples <- function(samples_filepath) {

  if (!file.exists(samples_filepath)) {
    stop(paste0("unable to find samples file ", samples_filepath))
  }
  
  samples <- read.csv(samples_filepath,  stringsAsFactors = FALSE)
  
  if (colnames(samples)[1] != "sample_id") {
    stop("first column of samples.csv should be sample_id")
  }
  
  if (any(duplicated(samples$sample_id))) {
    stop("sample_id should be unique")
  }
  
  if (ncol(samples) == 1) {
    stop("samples file includes no base FCS file designations")
  }
  
  samples$sample_id <- as.character(samples$sample_id)
  
  samples
}


#' Import FCS file
#'
#' Load data from an FCS file and convert it to a data frame with column names
#' set to a concatenation of fluorophore and marker.
#'
#' @param filepath Filepath of an FCS file. 
#' @return A data frame with expression data from an FCS file. If return_original is
#' true, the function will return the original flowFrame object in addition to the
#' data frame.
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


#' Generate population assignments for manually gated data
#'
#' Generate a list with keys representing sample IDs (ex: sample_1, sample_2).
#' Values are data frames with columns representing cell populations (ex: 
#' sample_1.basophil, sample_1.b_cell). Data frame cell values are booleans. 
#'
#' @param clustering_channels Vector containing clustering channels (ex: c("Er167Di_CCR7",
#' "Er168Di_CD3"))
#' @param manual_labeling_filepath Path to csv file with columns representing 
#' sample_id, base, and manually gated populations. Each row contains
#' a sample_id as well as the corresponding base FCS filename and FCS files for
#' manually gating populations. Note: The base FCS file contains data for all 
#' cells from the sample, regardless of manual gating population designation
#' @param data_dir Path to directory containing all FCS files.
#' @inheritParams .read_samples
#' @return A list illustrating the relationships between specific samples, cells,
#' and population assignments designated via manual gating. 
#' @export

generate_population_assignments <- function(clustering_channels, manual_labeling_filepath, samples_filepath, data_dir) {
  manual_labeling <- read.csv(manual_labeling_filepath, stringsAsFactors = FALSE)
  labels <- setdiff(names(manual_labeling), c("sample_id", "base"))
  single_sample_analysis_labels <- list()
  samples <- .read_samples(samples_filename)
  sample_ids <- samples$sample_id

  # Match manual gating labels to each sample.
  single_sample_labels <- lapply(sample_ids, function(sample_id) {
    cat(paste0(sample_id, "\n"))
    
    sample_labeling <- manual_labeling[manual_labeling$sample_id == sample_id, ]
    sample_base <- sample_labeling$base
    
    base_fcs_data <- .fcs_import_file(file.path(data_dir, sample_base))

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


#' Greedy search for optimal Average Overlap Frequency (AOF) for one channel.
#'
#' Given a cytometry data matrix and its clustering scheme, use a greedy search
#' strategy to find a partition of clusters that minimizes the AOF for each
#' marker.
#'
#' @param y_indices A list of indices vectors, one for each cluster.
#' @return Optimal AOF for x.
#' @export
#' @inheritParams calculateAof
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

#' Greedy search for optimal Average Overlap Frequency (AOF) values.
#'
#' Given a cytometry data matrix and its clustering scheme, use a greedy search
#' strategy to find a partition of clusters that minimizes the AOF for each
#' marker.
#'
#' @param fcs_data A numerical matrix (NxD) of acquired cytometry data. Each row
#' corresponds to a cell, each column to a channel.
#' @param y A factor vector of length N which includes the cluster assignment of
#' each cell.
#' @param channel_names A character vector which lists the columns (channels) in
#' fcs_data for which to calculate the AOF.
#' @param verbose Boolean. Should the function report progress to terminal.
#' @return A data frame with AOF values for each channel.
#' @export
#' @inheritParams calculateAof
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
