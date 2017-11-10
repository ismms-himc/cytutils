#' FCS channel renaming.
#'
#' A two-step utility for renaming FCS channel names and descriptions. This is
#' useful for situations where channels were misnamed across acquisitions (for
#' example, naming a channel CD8 in one sample and CD8a in another). See the
#' online documentation at TODO for instructions.
#'
#' @param path Directory where FCS files are located.
#' @param verbose If TRUE, script will report progress to console.
#' @export
channelRename <- function(path, verbose = TRUE) {
  filenames <- file.path(path, dir(path, "\\.fcs$"))
  if (length(filenames) == 0) stop("could not find any FCS files at given path")

  csv_filename <- file.path(path, "channel_rename.csv")
  if (file.exists(csv_filename)) {
    if (verbose) {
      message("Importing channel_rename.csv and exporting new FCS files")
    }

    channels <- read.csv(csv_filename, stringsAsFactors = FALSE)
    dest_path <- file.path(path, "channel_rename")
    if (!file.exists(dest_path)) dir.create(dest_path)
    renameFcsFileChannels(dest_path, filenames, channels, verbose = verbose)

    if (verbose) message("Export done")
  } else {
    if (verbose) message("Generating a new channel_rename.csv file")

    channels <- importChannelNames(filenames)

    # Export CSV.
    channels <- channels[order(channels$mass), ]
    channels$new_name <- channels$name
    channels$new_desc <- channels$desc
    write.csv(channels, csv_filename, row.names = FALSE)

    if (verbose) {
      message(paste0(
        "channel_rename.csv created. Update the file and re-run channelRename ",
        "to export new FCS files"
      ))
    }
  }
}

#' Import channel names and descriptions.
#'
#' Given a list of FCS file names, import channel names and descriptions across
#' all files.
#'
#' @return Data frame with unique channel names and descriptions, channel mass,
#' and whether that mass repeats across files with different name or
#' description.
#' @param filenames List of FCS file names.
#' @inheritParams channelRename
#' @export
importChannelNames <- function(filenames, verbose = TRUE) {
  # Read channel name + desc from FCS headers (don't load entire file) and
  # convert to a data frame.
  channels <- lapply(filenames, function(filename) {
    if (verbose) message(paste0("\t", filename))

    header <- flowCore::read.FCSheader(filename)
    header <- header[[1]]
    n_channels <- as.numeric(header[["$PAR"]])
    name <- unname(header[paste("$P", seq(n_channels), "N", sep = "")])
    desc <- unname(header[paste("$P", seq(n_channels), "S", sep = "")])

    data.frame(
      name = name,
      desc = desc,
      stringsAsFactors = FALSE
    )
  })
  channels <- do.call(rbind, channels)
  channels <- unique(channels)

  # Get mass from channel names.
  mass <- gsub("[^0-9]", "", channels$name)
  channels$mass <- as.integer(mass)
  channels <- channels[!is.na(channels$mass), ]
  channels <- channels[, c("mass", "name", "desc")]
  # Tag masses that have different names or descriptions across files.
  channels$dup <- channels$mass %in% channels$mass[duplicated(channels$mass)]

  channels
}

#' Rename channels in set of FCS files.
#'
#' Given a list of FCS files and a channel rename data frame, import each of the
#' files, rename their channel names and descriptions, and export them as a new
#' file.
#'
#' @param dest_path Directory where new FCS files are exported.
#' @param channels Channel rename data frame.
#' @param suffix New FCS file suffix.
#' @inheritParams importChannelNames
renameFcsFileChannels <- function(dest_path,
                                  filenames,
                                  channels,
                                  suffix = "renamed",
                                  verbose = TRUE) {
  # Only keep channels that require changes.
  channels <-
    channels[channels$name != channels$new_name |
               channels$desc != channels$new_desc, ]

  # Fix each file in turn.
  for (filename in filenames) {
    if (verbose) message(paste0("\t", filename))

    fcs <- flowCore::read.FCS(filename)
    # Merge file parameters with channel renames.
    params <- flowCore::pData(flowCore::parameters(fcs))
    params$keyword <- rownames(params)
    params <- merge(params, channels)

    # Rename each parameter in turn.
    desc <- description(fcs)
    for (row_idx in seq(nrow(params))) {
      row <- params[row_idx, ]
      # flowCore::write.FCS takes name from flowFrame column names.
      colnames(fcs)[colnames(fcs) == row$name] <- row$new_name
      # flowCore::write.FCS takes desc from flowFrame description field.
      desc[paste0(row$keyword, "S")] <- row$new_desc
    }

    # Update FlowFrame and export.
    description(fcs) <- desc
    filename <- basename(filename)
    if (suffix != "") {
      dest_filename <-
        file.path(dest_path, paste0(filename, ".", suffix, ".fcs"))
    } else {
      dest_filename <-
        file.path(dest_path, filename)
    }
    if (verbose) message("\t--> ", dest_filename)
    flowCore::write.FCS(fcs, dest_filename)
  }
}
