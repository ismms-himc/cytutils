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
  exprs <-
    apply(exprs, 2, function(x) (x - min(x)) / (quantile(x, 0.99) - min(x)))
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
    Label = labels,
    stringsAsFactors = FALSE
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
                                    key,
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
  diffs <- unlist(lapply(seq(n_pos - 1), function(pos) {
    exprs_list$exprs_sorted[, pos] - exprs_list$exprs_sorted[, pos + 1]
  }))
  bc_separation_dist <-
    with(exprs_list, exprs_sorted[, n_pos] - exprs_sorted[, n_pos + 1])
  diff_thresh <- quantile(diffs, threshold_per)
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

#' Export debarcoded FCS files.
#'
#' Given an FCS file, export a separate FCS file for each of the event labels
#' (including unlabeled). If any additional parameters were used in the labeling
#' (such as Mahalanobis Ratio) they will be included as columns in the FCS data.
#'
#' @param path_prefix Prefix to path where files will be exported. Each file will
#' be named [prefix_path].[label].fcs.
#' @inheritParams debarcoderPrepareFcs
#' @inheritParams debarcoderUnlabelEvents
#' @export
debarcoderExportDebarcodedFcs <- function(path_prefix, fcs, labels) {
  cols <- setdiff(colnames(labels), "Label")
  codes <- unique(labels$Label)

  # Export each code in turn.
  for (code in codes) {
    code_indices <- which(labels$Label == code)
    code_fcs <- fcs[code_indices, ]

    # Add any columsn that were involved in labeling.
    if (length(cols) > 0) {
      e <- flowCore::exprs(code_fcs)

      for (col in cols) {
        e <- cbind(e, labels[code_indices, col])
        colnames(e)[ncol(e)] <- col
      }

      code_fcs <- flowCore::flowFrame(e, description = code_fcs@description)
    }

    flowCore::write.FCS(code_fcs, paste0(path_prefix, ".", code, ".fcs"))
  }
}

#' Debarcoder diagnostic plots.
#'
#' Given an FCS file, generate a separate scatter plot of Mahalanobis ratio
#' versus barcoding separation distance.
#'
#' @param path_prefix Either NULL (do not export figures) or prefix to path
#' where files will be exported. Each file will be named
#' [prefix_path].[label].jpg.
#' @inheritParams debarcoderPrepareFcs
#' @inheritParams debarcoderUnlabelEvents
#' @importFrom ggplot2 ggplot aes geom_point scale_x_log10 scale_y_continuous labs theme ggsave
#' @return A list of ggplot objects corresponding to the figures.
#' @export
debarcoderPlots <- function(path_prefix, fcs, labels) {
  if (!("BcSepDist" %in% colnames(labels))) stop("labels missing BcSepDist")
  if (!("MahalRatio" %in% colnames(labels))) stop("labels missing MahalRatio")

  codes <- unique(labels$Label)

  xlim <- c(1, max(labels$MahalRatio))

  figs <- lapply(codes, function(code) {
    code_indices <- which(labels$Label == code)
    ggplot(labels[code_indices, ], aes(x = MahalRatio, y = BcSepDist)) +
      geom_point(size = 0, alpha = 0.2) +
      scale_x_log10(limits = xlim) +
      scale_y_continuous(limits = c(0, 1)) +
      labs(title = paste0("Debarcoding diagnostic for ", code),
           x = "log10(Mahalanobis Ratio)",
           y = "Barcoding Separation Distance") +
      theme(aspect.ratio = 1)
  })
  names(figs) <- codes

  # Export figures as JPGs if path_prefix exists.
  if (!is.null(path_prefix)) {
    for (code in codes) {
      ggsave(paste0(path_prefix, ".", code, ".jpg"),
             figs[[code]], width = 4, height = 4)
    }
  }

  figs
}

#' Debarcode an FCS file.
#'
#' Import FCS file and barcoding key and run it through the debarcoding
#' workflow: label events using Zunder et al., unlabel events according to
#' barcoding separation distance and Mahalanobis ratio heuristics, and export
#' each label to a separate file along with an accompanying diagnostic plot.
#'
#' Please note that this function ends with a series of warnings due to
#' flowCore's write.FCS.
#'
#' @param fcs_file_path FCS file path.
#' @param key_file_path CSV barcoding key file path.
#' @param export_figures If TRUE, the function will export diagnostic plots.
#' @param verbose If TRUE, the function will message the console with updates.
#' @export
debarcode <- function(fcs_file_path,
                      key_file_path,
                      export_figures = TRUE,
                      verbose = TRUE) {
  if (!file.exists(fcs_file_path)) {
    stop(paste0("unable to find ", fcs_file_path))
  }
  if (!file.exists(key_file_path)) {
    stop(paste0("unable to find ", key_file_path))
  }

  key <- cytutils::debarcoderImportKey(key_file_path)
  if (verbose) message("importing FCS file ...")
  fcs <- flowCore::read.FCS(fcs_file_path)
  if (verbose) message("debarcoding ...")
  exprs_list <- cytutils::debarcoderPrepareFcs(fcs, key)
  labels <- cytutils::debarcoderLabelEvents(exprs_list, key)
  if (verbose) message("unlabeling events using heuristics ...")
  labels <- cytutils::debarcoderUnlabelEvents(exprs_list, labels, key)
  if (verbose) message("exporting debarcoded FCS files ...")
  cytutils::debarcoderExportDebarcodedFcs(fcs_file_path, fcs, labels)
  if (export_figures) {
    if (verbose) message("exporting plots ...")
    cytutils::debarcoderPlots(fcs_file_path, fcs, labels)
  }

  return()
}
