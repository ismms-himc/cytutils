


#' Calculate the Average Overlap Frequency (AOF).
#'
#' Calculate the Average Overlap Frequency (AOF) for a cytometry data matrix.
#'
#' @param fcs_data A numerical matrix (NxD) of acquired cytometry data. Each row
#' corresponds to a cell, each column to a channel.
#' @param y A factor vector of length N with cell assignments.
#' @param channel_names A character vector which lists the columns (channels) in
#' fcs_data for which to calculate the AOF.
#' @param cofactor A numeric. If specified, fcs_data will be transformed using
#' inverse hyperbolic sin (arcsinh) with this cofactor.
#' @return A data frame with AOF values for each channel.
#' @export
calculateAof <- function(fcs_data,
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

}
