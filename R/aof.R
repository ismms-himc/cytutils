#' Average Overlap Frequency (AOF).
#'
#' Calculate the Average Overlap Frequency (AOF) for a single cytometry channel.

#' @param x Numeric vector corresponding to one channel in cytometry data.
#' @param pos_indices Indices of cells positive for this channel.
#' @param width Width of high threshold of negative populationa and low
#' threshold of positive population.
#' @param cofactor If supplied, x will be transformed using inverse hyperbolic
#' sin with given cofactor.
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

#' TODO
.greedyChannelAof <- function(x, y, cofactor = NULL) {
  n <- length(x)
  df <- tibble::tibble(x = x, y = NA)
  # TODO fix this
  for (class_label in names(class_indices)) {
    df$class_label[class_indices[[class_label]]] <- class_label
  }
  df_stats <- df %>%
    dplyr::group_by(class_label) %>%
    dplyr::summarize(freq = n() / n, mean_v = mean(v)) %>%
    dplyr::arrange(mean_v)
  df_stats$cum_sum_freq = cumsum(df_stats$freq)

  find_x_y <- function(df_stats, cum_sum_thresh) {
    x_class_labels <- as.character(
      dplyr::filter(df_stats, cum_sum_freq <= cum_sum_thresh)$class_label)
    x <- df$v[unlist(class_indices[x_class_labels])]
    y_class_labels <-
      as.character(setdiff(df_stats$class_label, x_class_labels))
    y <- df$v[unlist(class_indices[y_class_labels])]

    list(x = x, y = y)
  }

  compute_aof <- "aof" %in% metrics
  compute_si <- "si" %in% metrics

  # cum sum thresholds for search
  search_thresh <- c(
    seq(0.2, 0.5, 0.05),
    seq(0.5, 0.8, 0.03),
    seq(0.8, 0.98, 0.02), 0.99
  )

  qc_metrics <- lapply(search_thresh, function(cum_sum_thresh) {
    xy <- find_x_y(df_stats, cum_sum_thresh)
    x <- xy$x
    y <- xy$y

    cum_sum_thresh_aof <- Inf
    cum_sum_thresh_si <- -Inf

    if (length(x) > 0 & length(y) > 0) {
      if (compute_aof)  cum_sum_thresh_aof <- aof(x, y)
      if (compute_si) {
        cum_sum_thresh_si <- si(x, y)
        if (is.nan(cum_sum_thresh_si)) cum_sum_thresh_si <- 0
      }
    }

    dplyr::data_frame(
      cum_sum_thresh = cum_sum_thresh,
      aof = cum_sum_thresh_aof,
      si = cum_sum_thresh_si
    )
  }) %>% dplyr::bind_rows()

  # return minimum aof and maximum si
  output <- list()

  if (compute_aof) {
    min_aof_index <- min(which(qc_metrics$aof == min(qc_metrics$aof)))
    min_aof_xy <- find_x_y(df_stats, qc_metrics$cum_sum_thresh[[min_aof_index]])
    output$aof <- list(
      value = qc_metrics$aof[min_aof_index],
      x = min_aof_xy$x,
      y = min_aof_xy$y
    )
  }

  if (compute_si) {
    max_si_index <- min(which(qc_metrics$si == max(qc_metrics$si)))
    max_si_xy <- find_x_y(df_stats, qc_metrics$cum_sum_thresh[[max_si_index]])
    output$si <- list(
      value = qc_metrics$si[max_si_index],
      x = max_si_xy$x,
      y = max_si_xy$y
    )
  }

  output
}


#' Greedy search for optimal Average Overlap Frequency (AOF) values.
#'
#' Given a cytometry data matrix and its clustering scheme, use a greedy search
#' strategy to find a partition of clusters that minimizes the AOF for each
#' marker.
#'
#' @param fcs_data A numerical matrix (NxD) of acquired cytometry data. Each row
#' corresponds to a cell, each column to a channel.
#' @param y A factor vector of length N which includes the cluster of each cell.
#' @param channel_names A character vector which lists the columns (channels) in
#' fcs_data for which to calculate the AOF.
#' @param cofactor A numeric. If specified, fcs_data will be transformed using
#' inverse hyperbolic sin (arcsinh) with this cofactor.
#' @param verbose Boolean. Should the function report progress to terminal.
#' @return A data frame with AOF values for each channel.
#' @export
greedyCytometryAof <- function(fcs_data,
                               y,
                               channel_names = colnames(fcs_data),
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

  # Run greedy algorithm over each channel.
  aof_values <- lapply(channel_names, function(channel_name) {
    if (verbose) message(channel_name)
    x <- fcs_data[[channel_name]]
    data.frame(
      ChannelName = channel_name,
      Aof = .greedyChannelAof(x, y, cofactor)
    )
  })

  do.call(rbind, aof_values)
}
