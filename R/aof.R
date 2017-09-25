#' Average Overlap Frequency (AOF).
#'
#' Calculate the Average Overlap Frequency (AOF) for a single cytometry channel.
#'
#' @param x Numeric vector corresponding to one channel in cytometry data.
#' @param pos_indices Indices of cells positive for this channel.
#' @param width Width of high threshold of negative populationa and low
#' threshold of positive population.
#' @param cofactor If supplied, data will be transformed using inverse
#' hyperbolic sin with given cofactor.
#' @return The AOF between positive and negative populations for x.
#' @export
calculateAof <- function(x, pos_indices, width = 0.05, cofactor = NULL) {
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
  neg_indices <- setdiff(seq_along(x), pos_indices)
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
    }

    calculateAof(x, pos_indices, width = width, cofactor = cofactor)
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
                               verbose = FALSE) {
  if (!is.factor(y)) {
    stop("y should be a factor vector.")
  }

  if (length(y) != nrow(fcs_data)) {
    stop("y length should be equal to number of rows in fcs_data")
  }

  if (!all(channel_names %in% colnames(fcs_data))) {
    stop("channel_names should all be columns of fcs_data")
  }

  if (any(is.nan(fca_data[, channel_names]))) {
    stop("cannot calculate AOF when data includes nan values")
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
    x <- fcs_data[[channel_name]]
    data.frame(
      ChannelName = channel_name,
      Aof = .greedyChannelAof(x, y_indices, width = width, cofactor = cofactor)
    )
  })

  do.call(rbind, aof_values)
}
