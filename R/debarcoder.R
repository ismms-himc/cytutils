.verifyExprsList <- function(exprs_list) {
  if (!is.list(exprs_list)) stop("exprs_list is not a list")
  if (is.null(exprs_list$exprs)) stop("exprs_list is missing exprs field")
  if (is.null(exprs_list$exprs_sorted)) {
    stop("exprs_list is missing exprs_sorted field")
  }
}

.verifyKey <- function(key) {
  if (!is.list(key)) stop("key is not a list")
  if (is.null(key$key)) stop("key is missing key field")
  if (is.null(key$channels)) stop("key is missing channels field")
  if (is.null(key$n_pos_channels)) stop("key is missing n_pos_channels field")
}

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
#' Transform expression data using cofactor, scale to [0..1] range, match to
#' barcoding channels, and sort each event separately.
#'
#' @param fcs FCS flow frame.
#' @param key Barcoding key data frame.
#' @param cofactor Before scaling, data will be transformed using inverse
#' hyperbolic sin with given cofactor.
#' @return Expression data after transformation, scaling, and matching, and an
#' additional copy after sorting.
#' @export
debarcoderPrepareFcs <- function(fcs, key, cofactor = 10) {
  if (class(fcs) != "flowFrame") stop("fcs is not a flowFrame object")
  .verifyKey(key)

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
  exprs <- exprs[, key$channels]
  exprs_sorted <- t(apply(exprs, 1, sort, decreasing = TRUE))

  list(
    exprs = exprs,
    exprs_sorted = exprs_sorted
  )
}


#' Label events.
#'
#' For each event, label it with the barcoding code that matches that event's
#' top N highest intensity channels.
#'
#' @param exprs_list A list with expression and sorted expression matrices.
#' @param unlabeled_label Label for events that don't match any key codes.
#' @inheritParams debarcoderPrepareFcs
#' @return A data frame with a row for every event and a single Label column.
#' @export
debarcoderLabelEvents <- function(exprs_list,
                                  key,
                                  unlabeled_label = "unlabeled") {
  .verifyExprsList(exprs_list)
  .verifyKey(key)

  event_thresh <- exprs_list$exprs_sorted[, key$n_pos_channels]
  exprs_bool <- as.data.frame(exprs_list$exprs >= event_thresh)
  exprs_bool$Index <- seq(nrow(exprs_bool))
  exprs_bool <- merge(exprs_bool, key$key, all.x = TRUE)
  labels <- exprs_bool$code[order(exprs_bool$Index)]
  labels[is.na(labels)] <- unlabeled_label

  data.frame(
    Label = labels
  )
}

#' Unlabel events.
#'
#' Identify a set of suspect events whose barcoding separation distance is below
#' a given threshold. Then, unlabel events if their Mahalanobis ratio is above
#' the threshold set by the suspect events.
#'
#' Please consult the package documentation for a complete explanation of the
#' heuristics involved (TODO).
#'
#' @param labels Data frame with an event label column.
#' @param threshold_per Percentile for setting distance and ratio thresholds.
#' @inheritParams debarcoderLabelEvents
#' @return A data frame with a row for every event. Columns are label, barcoding
#' separation distance, and Mahalanobis ratio.
#' @export
debarcoderUnlabelEvents <- function(exprs_list,
                                    labels,
                                    threshold_per = 0.99,
                                    unlabeled_label = "unlabeled") {
  .verifyExprsList(exprs_list)
  if (!is.data.frame(labels)) stop("labels is not a data frame")
  if (is.null(labels$Label)) stop("labels is missing the Label column")
  if (!is.null(labels$BcSepDist) || !is.null(labels$MahalRatio)) {
    stop("labels has already been unlabeled")
  }

  # Identify suspect events.
  n_pos <- key$n_pos_channels
  diff_1_2 <- with(exprs_list, exprs_sorted[, 1] - exprs_sorted[, 2])
  bc_separation_dist <-
    with(exprs_list, exprs_sorted[, n_pos] - exprs_sorted[, n_pos + 1])
  diff_thresh <- quantile(diff_1_2, threshold_per)
  suspect_indices <- which(bc_separation_dist < diff_thresh)

  # Calculate Mahalanobis ratio.
  codes <- sort(setdiff(labels$Label, unlabeled_label))

  indices <- lapply(codes, function(code) which(labels$Label == code))
  names(indices) <- codes
  centers <- lapply(indices, function(idx) colMeans(exprs_list$exprs[idx, ]))
  covs <- lapply(indices, function(idx) cov(exprs_list$exprs[idx, ]))

  dists <- lapply(codes, function(code) {
    mahalanobis(exprs_list$exprs, centers[[code]], covs[[code]])
  })
  dists <- do.call(cbind, dists)
  dists_sorted <- t(apply(dists, 1, sort))
  mahal_ratio <- dists_sorted[, 2] / dists_sorted[, 1]

  # Unlabel events whose Mahalanobis ratio is in the suspect range.
  ratio_thresh <- quantile(mahal_ratio[suspect_indices], threshold_per)
  unlabel_indices <- which(mahal_ratio <= ratio_thresh)
  labels$Label[unlabel_indices] <- unlabeled_label

  labels$BcSepDist <- bc_separation_dist
  labels$MahalRatio <- mahal_ratio
  labels
}
