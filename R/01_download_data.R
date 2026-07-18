# Download GEO supplementary files and list what was actually retrieved.
# Do not assume file format here — inspect the listing before writing loaders.

source(file.path(
  if (basename(getwd()) == "R") ".." else ".",
  "R", "00_config.R"
))

if (!requireNamespace("GEOquery", quietly = TRUE)) {
  stop(
    "Package 'GEOquery' is required.\n",
    "Install with:\n",
    "  if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager')\n",
    "  BiocManager::install('GEOquery')"
  )
}

suppressPackageStartupMessages(library(GEOquery))

#' Recursively list files under a directory with sizes.
list_supp_tree <- function(root_dir) {
  if (!dir.exists(root_dir)) {
    return(data.frame(
      path = character(),
      bytes = numeric(),
      stringsAsFactors = FALSE
    ))
  }
  paths <- list.files(root_dir, recursive = TRUE, full.names = TRUE, include.dirs = FALSE)
  if (length(paths) == 0) {
    return(data.frame(path = character(), bytes = numeric(), stringsAsFactors = FALSE))
  }
  info <- file.info(paths)
  data.frame(
    path = normalizePath(paths, winslash = "/", mustWork = FALSE),
    bytes = as.numeric(info$size),
    stringsAsFactors = FALSE
  )
}

#' Expand common archive types downloaded by GEOquery (in place).
expand_archives <- function(root_dir) {
  archives <- list.files(
    root_dir,
    pattern = "\\.(tar\\.gz|tgz|tar|zip)$",
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = TRUE
  )
  for (arc in archives) {
    message("Expanding archive: ", arc)
    out_dir <- dirname(arc)
    tryCatch(
      {
        if (grepl("\\.zip$", arc, ignore.case = TRUE)) {
          utils::unzip(arc, exdir = out_dir)
        } else {
          utils::untar(arc, exdir = out_dir)
        }
      },
      error = function(e) {
        warning("Failed to expand ", arc, ": ", conditionMessage(e))
      }
    )
  }
}

download_one_gse <- function(gse) {
  dest <- raw_gse_dir(gse)
  if (!dir.exists(dest)) {
    dir.create(dest, recursive = TRUE, showWarnings = FALSE)
  }

  message("\n========================================")
  message("Downloading supplementary files for ", gse)
  message("Destination: ", dest)
  message("========================================")

  result <- tryCatch(
    {
      GEOquery::getGEOSuppFiles(
        GEO = gse,
        makeDirectory = FALSE,
        baseDir = dest,
        fetch_files = TRUE
      )
    },
    error = function(e) {
      message("ERROR: getGEOSuppFiles() failed for ", gse, ": ", conditionMessage(e))
      NULL
    }
  )

  if (is.null(result)) {
    message("No files retrieved for ", gse, " (download failed).")
    return(invisible(NULL))
  }

  # Expand archives so 02_preprocess can see mtx/h5/tsv contents
  expand_archives(dest)

  listing <- list_supp_tree(dest)
  listing_path <- file.path(dest, "file_listing.txt")
  write.table(
    listing,
    file = listing_path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  message("\n--- File listing for ", gse, " (", nrow(listing), " files) ---")
  if (nrow(listing) == 0) {
    message("(empty — check GEO accession / network)")
  } else {
    sizes_mb <- sprintf("%.2f MB", listing$bytes / 1e6)
    for (i in seq_len(nrow(listing))) {
      message(sprintf("  [%s] %s", sizes_mb[i], listing$path[i]))
    }
  }
  message("Listing also written to: ", listing_path)
  message(
    "Inspect this listing before relying on 02_preprocess.R loaders ",
    "(10x mtx / .h5 / flat matrix formats vary by GSE)."
  )

  invisible(listing)
}

main <- function() {
  message("Stage 01: download GEO supplementary files")
  for (ds in DATASETS) {
    download_one_gse(ds$gse)
  }
  message("\nStage 01 complete. Review the listings above, then run 02_preprocess.R.")
}

main()
