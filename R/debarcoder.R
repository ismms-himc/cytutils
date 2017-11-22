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

#' Sort barcoding channel intensities.
#'
#' Given an FCS flow frame, get intensity (expression) values for the barcoding
#' channels and sort them for each event.
#'
#' @param fcs FCS flow frame.
#' @param key Barcoding key data frame.
#' @return Matrix of intensity values, sorted by event.
debarcoderSortExprs <- function(fcs, key) {
  # Match FCS channels to barcoding channels.
  missing_channels <- setdiff(key$channels, fcs@parameters@data$name)
  if (length(missing_channels) > 0) {
    stop("FCS data is missing the following barcoding channels: ",
         paste(missing_channels, collapse = ", "))
  }

  exprs <- flowCore::exprs(fcs)[, key$channels]
  t(apply(exprs, 1, sort, decreasing = TRUE))
}
