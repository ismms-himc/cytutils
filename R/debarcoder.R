#' Import a barcoding key.
#'
#' Import a barcoding key from a CSV file. The CSV header should be of the
#' format:
#'
#' code,channel_name_1,channel_name_2,...
#'
#' Channel names should be explicit (such as Pd102Di). Channel values should be
#' denoted as 1 (channel positive for code) or 0 (channel negative for code).
#'
#' @param filename CSV file from which to import the key.
#' @return List with barcoding key, barcoding channel names, and number of
#' positive channels for each code.
#' @export
debarcoderImportKey <- function(filename) {
  key <- read.csv(filename, stringsAsFactors = FALSE)

  n <- nrow(key)
  channels <- colnames(key)[2:ncol(key)]
  channel_data <- key[, channels]
  colnames(key)[1] <- tolower(colnames(key)[1])

  if (colnames(key)[1] != "code") {
    stop("first column of key should be named \"code\"")
  }

  if (!all(channel_data == 0 | channel_data == 1)) {
    stop("channel values are not 0 or 1")
  }

  if (nrow(unique(channel_data)) != n) {
    stop("barcoding key cannot include duplicate code channel combinations")
  }

  if (length(unique(key$code)) != n) {
    stop("barcoding key cannot include duplicate code names")
  }

  n_pos_channels <- unique(rowSums(channel_data))
  if (length(n_pos_channels) > 1) {
    stop("codes have a variable number of positive channels")
  }

  # Convert channel values to boolean.
  key[, channels] <- channel_data == 1

  list(
    key = key,
    channels = channels,
    n_pos_channels = n_pos_channels
  )
}

#' Prepare FCS data for debarcoding.
#'
#' Transform expression data using cofactor, scale to [0..1] range, and match to
#' barcoding channels.
#'
#' @param fcs FCS flow frame.
#' @param key Barcoding key data frame.
#' @return Expression data after transformation, scaling, and matching.
debarcoderPrepareFcs <- function(fcs, key, cofactor = 10) {
  # Verify FCS has all key channels.
  missing_channels <- setdiff(key$channels, fcs@parameters@data$name)
  if (length(missing_channels) > 0) {
    stop("FCS data is missing the following barcoding channels: ",
         paste(missing_channels, collapse = ", "))
  }

  # Transform and scale.
  exprs <- flowCore::exprs(fcs)
  exprs <- asinh(exprs / cofactor)
  exprs <- apply(exprs, 2, function(x) (x - min(x)) / (max(x) - min(x)))
  # Match FCS channels to barcoding channels.
  exprs[, key$channels]
}
