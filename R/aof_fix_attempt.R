# sum(pos <= neg_high, na.rm = TRUE)
# TODO: added default neg_indices to account for when greedyAOF functions call this function
test_calculateAof <- function(x, pos_indices, neg_indices = setdiff(seq_along(x), pos_indices), width = 0.05, cofactor = 5) {
  if (length(pos_indices) == 0) {
    stop("no cells in positive population")
  }

  if (length(neg_indices) == 0) {
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
    # TODO: switched neg to pos (BAD IDEA)
    # neg_high <- qnorm(1 - width, mean(pos), sd(pos)) # VERY BAD AOF
  }

  if (length(pos) == 1) {
    pos_low <- pos
  } else {

    # pos_low <- quantile(pos, probs = c(.05))
    # quantile_pos <- quantile(pos, probs = seq(0,1,.05))
    # pos_low <- quantile_pos['5%']

    # pos_cleaned <- pos[!pos %in% boxplot(pos)$out]
    # Note: qnorm is the invere of pnorm
    pos_low <- qnorm(width, mean(pos), sd(pos))

    # TODO: switched pos to neg (BAD IDEA)
    # pos_low <- qnorm(width, mean(pos), sd(pos)) # VERY BAD AOF
    # TODO: added abs(width -1) below
    # pos_low <- qnorm(abs(width - 1), mean(neg), sd(neg)) AOF = .47
  }

  # Calculate AOF.
  # mean(c(mean(pos <= neg_high), mean(neg >= pos_low)))
  
  num_pos_below <- sum(pos <= neg_high, na.rm = TRUE)
  num_neg_above <- sum(neg >= pos_low, na.rm = TRUE)
  # length(neg >= pos_low)
  num_pos <- length(pos)
  num_neg <- length(neg)

  0.5 * ((num_neg_above / num_neg) + (num_pos_below / num_pos))
}

test_calculateAof(x, t_cell_indices, non_t_cell_indices) # =>  0.002309506

# `pos` is a vector of Er168Di intensity values (positive floats) for t cells
# `length(pos)` #=> 154605
# `mean(pos)` #=> 385.6191
# `sd(pos)` #=> 250.2473
# `width <- 0.05`
# `qnorm(width, mean(pos), sd(pos))` #=> -26.00106


widths <- seq(0, 1, 0.01)
qnorms <- sapply(widths, function(width) {
  qnorm(width, mean(pos), sd(pos))
})

pnorms <- sapply(widths, function(width) {
  pnorm(width, mean(pos), sd(pos))
})



# widths = (0..1).step(0.01) # numbers from 0 to 1 with steps of 0.01 between values
# qnorms = widths.map { |width| qnorm(width, mean, sd) }
# plot(widths, qnorms)


# -----------------------------------------------------------


test_greedyCytometryAof <- function(fcs_data,
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
      Aof = .test_greedyChannelAof(x, y_indices, width = width, cofactor = cofactor)
    )
  })

  do.call(rbind, aof_values)
}


.test_greedyChannelAof <- function(x, y_indices, width = 0.05, cofactor = NULL) {
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
      aof <- test_calculateAof(x, pos_indices, width = width, cofactor = cofactor)
    }

    aof
  })

  # Return the minimal AOF found.
  min(unlist(aof_values))
}


test_greedyCytometryAof(base_fcs_data@exprs, cell_assignments_ordered, channel_names)

# ChannelName         Aof
# 1     Er168Di 0.130286586
# 2     Nd142Di 0.010353629
# 3     Gd158Di 0.006306912
# 4     Dy161Di 0.084751372