#' FCS channel renaming.
#'
#' A two-step utility for renaming FCS channel names and descriptions. This is
#' useful for situations where channels were misnamed across acquisitions (for
#' example, naming a channel CD8 in one sample and CD8a in another). See the
#' Readme at "FCS Channel Renaming" for instructions.
#'
#' @param path Directory where FCS files are located.
#' @param dup_handling How to handle masses with duplicates name or desc.
#' Options are "stop", "message", and NULL.
#' @param verbose If TRUE, script will report progress to console.
#' @param na.mass.rm Remove <NA> mass (TRUE by default)
#' @param ignore.mass Ignore the mass computation (FALSE by default); if TRUE
#' the mass column is assigned <NA>
#' @export
channelRename <- function(path, dup_handling = "message", verbose = TRUE, na.mass.rm = TRUE, ignore.mass = FALSE) {
  filenames <- file.path(path, dir(path, "\\.fcs$"))
  if (length(filenames) == 0) stop("could not find any FCS files at given path")

  if (!is.null(dup_handling)) {
    if (!(dup_handling %in% c("message", "stop"))) {
      stop("dup_handling not in allowed values (NULL, message, or stop)")
    }
  }

  csv_filename <- file.path(path, "channel_rename.csv")
  if (file.exists(csv_filename)) {
    if (verbose) {
      message("Importing channel_rename.csv and exporting new FCS files")
    }

    channels <- read.csv(csv_filename, stringsAsFactors = FALSE)
    channels$name[is.na(channels$name)] <- ""
    channels$desc[is.na(channels$desc)] <- ""

    if (!is.null(dup_handling)) {
      dups <- .findMassDups(channels)
      if (length(dups) > 0) {
        msg <-
          paste0("the following masses have duplicate new_name or new_desc: ",
                 paste(dups, collapse = ", "))
        if (dup_handling == "message") {
          message(msg)
        } else if (dup_handling == "stop") {
          stop(msg)
        }
      }
    }

    dest_path <- file.path(path, "channel_rename")
    if (!file.exists(dest_path)) dir.create(dest_path)
    renameFcsFileChannels(dest_path,
                          filenames,
                          channels,
                          verbose = verbose)

    if (verbose) message("Export done")
  } else {
    if (verbose) message("Generating a new channel_rename.csv file")

    channels <- importChannelNames(filenames, na.mass.rm, ignore.mass)

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
#' @param na.mass.rm Remove <NA> mass (TRUE by default)
#' @param ignore.mass Ignore the mass computation (FALSE by default); if TRUE
#' the mass column is assigned <NA>
#' @export
importChannelNames <- function(filenames, verbose = TRUE, na.mass.rm = TRUE, ignore.mass = FALSE) {
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
      filename = filename,
      stringsAsFactors = FALSE
    )
  })
  channels <- do.call(rbind, channels)
  channels <- channels[, c("name", "desc")]
  channels <- unique(channels)

  # Get mass from channel names.
  mass <- gsub("[^0-9]", "", channels$name)
  # If mass is ignored, update
  if (ignore.mass) {
    mass <- NA
    na.mass.rm = FALSE
  }
  channels$mass <- as.integer(mass)
  # NA removal is optionnal but TRUE by default
  if (na.mass.rm)
    channels <- channels[!is.na(channels$mass), ]
  channels <- channels[, c("mass", "name", "desc")]
  # Tag masses that have different names or descriptions across files.
  channels$dup <- channels$mass %in% channels$mass[duplicated(channels$mass)]
  # If mass is ignored, name is looked for duplication
  if (ignore.mass)
    channels$dup <- channels$name %in% channels$name[duplicated(channels$name)]

  channels
}

.findMassDups <- function(channels) {
  # Find masses with duplicate descriptions or names.
  desc_dups <- unique(channels[, c("mass", "new_desc")])
  name_dups <- unique(channels[, c("mass", "new_name")])

  c(desc_dups[duplicated(desc_dups$mass), "mass"],
    name_dups[duplicated(name_dups$mass), "mass"])
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
    params <- fcs@parameters@data
    params$name[is.na(params$name)] <- ""
    params$desc[is.na(params$desc)] <- ""

    params$keyword <- rownames(params)
    params <- merge(params, channels)

    filename <- basename(filename)
    if (suffix != "") {
      dest_filename <-
        file.path(dest_path, paste0(filename, ".", suffix, ".fcs"))
    } else {
      dest_filename <- file.path(dest_path, filename)
    }

    # This if clause prevents bug from occuring if the current fcs file already 
    # has appropriate channel names and descriptions.
    if (nrow(params) == 0) {
      if (verbose) message("\t--> ", dest_filename)
      flowCore::write.FCS(fcs, dest_filename)
    } else {
      # Rename each parameter in turn.
      desc <- fcs@description
      for (row_idx in seq(nrow(params))) {
        row <- params[row_idx, ]
        # Cover all possible venues for name, since it's unclear how flowCore
        # works. The github source says it's colnames(fcs@exprs) but sometimes
        # that does not work.
        colnames(fcs@exprs)[colnames(fcs@exprs) == row$name] <- row$new_name
        fcs@parameters@data[row$keyword, "name"] <- row$new_name
        fcs@parameters@data[row$keyword, "desc"] <- row$new_desc
        desc[paste0(row$keyword, "N")] <- row$new_name
        desc[paste0(row$keyword, "S")] <- row$new_desc
      }

      # Mark duplicated descriptions.
      descs <- fcs@parameters@data$desc
      descs <- descs[!is.na(descs)]
      dup_desc_keywords <- names(descs)[which(duplicated(descs))]
      for (keyword in dup_desc_keywords) {
        keyword <- gsub("S", "", keyword)
        new_desc <-
          paste0(fcs@parameters@data[keyword, "desc"], "_",
                 gsub("\\$", "", keyword))
        fcs@parameters@data[keyword, "desc"] <- new_desc
        desc[paste0(keyword, "S")] <- new_desc
      }

      # Update FlowFrame and export.
      fcs@description <- desc

      if (verbose) message("\t--> ", dest_filename)
      flowCore::write.FCS(fcs, dest_filename)
    }
  }
}