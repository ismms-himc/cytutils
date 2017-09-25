

#' Average Overlap Frequency (AOF).
#'
#' Calculate the Average Overlap Frequency (AOF) for a single cytometry channel.

#' @param x Numeric vector corresponding to one channel in cytometry data.
#' @param pos_indices Indices of cells positive for this channel.
#' @param width Width of high threshold of negative populationa and low
#' threshold of positive population.
#' @return The AOF between positive and negative populations for x.
#' @export
calculateAof <- function(x, pos_indices, width = 0.05) {
  if (max(pos_indices) > length(x)) {
    stop("positive population indices include values higher than length of x")
  }
  if (length(pos_indices) == 0) {
    stop("no cells in positive population")
  }
  if (length(pos_indices) == length(x)) {
    stop("no cells in negative population")
  }

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


#' Cytometry calculation of Average Overlap Frequency (AOF).
#'
#' Calculate the Average Overlap Frequency (AOF) for a cytometry data matrix.
#'
#' @param fcs_data A numerical matrix (NxD) of acquired cytometry data. Each row
#' corresponds to a cell, each column to a channel.
#' @param pos_indices A list of
#' @param channel_names A character vector which lists the columns (channels) in
#' fcs_data for which to calculate the AOF.
#' @param cofactor A numeric. If specified, fcs_data will be transformed using
#' inverse hyperbolic sin (arcsinh) with this cofactor.
#' @return A data frame with AOF values for each channel.
#' @export
calculateCytometryAof <- function(fcs_data,
                                  y,
                                  channel_names = colnames(fcs_data),
                                  cofactor = NULL) {
  if (!is.factor(y)) {
    stop("y should be a factor vector.")
  }

  if (length(y) != nrow(fcs_data)) {
    stop("y length should be equal to number of rows in fcs_data")
  }

  if (!all(channel_names %in% colnames(fcs_data))) {
    stop("channel_names should all be columns of fcs_data")
  }

  x <- fcs_data[, channel_names]

  if (!is.null(cofactor)) {
    x <- asinh(x / cofactor)
  }

  if (any(is.nan(x))) {
    stop("cannot calculate AOF when data includes nan values")
  }

}
