# Export scGFT stats + per-cluster synthesis counts as JSON for a separate frontend.
# Never invent accuracy/deviation — only serialize values produced by stage 03.

source(file.path(
  if (basename(getwd()) == "R") ".." else ".",
  "R", "00_config.R"
))

if (!requireNamespace("jsonlite", quietly = TRUE)) {
  stop(
    "Package 'jsonlite' is required.\n",
    "Install with: install.packages('jsonlite')"
  )
}

suppressPackageStartupMessages(library(jsonlite))

#' Build JSON array of {step, cellsSynthesized, accuracy, deviation}.
#'
#' Overall step carries real accuracy (statsScGFT) and deviation (parsed from
#' RunScGFT). Per-cluster steps carry synthesized cell counts; accuracy and
#' deviation are NULL there because statsScGFT returns a single global accuracy
#' and does not return per-cluster metrics.
build_export_rows <- function(result) {
  if (is.null(result$stats) || is.null(result$stats$accuracy)) {
    stop(
      "Missing stats$accuracy in results for ", result$gse,
      " — refusing to invent a value."
    )
  }

  accuracy <- as.numeric(result$stats$accuracy)
  if (length(accuracy) != 1 || is.na(accuracy)) {
    stop(
      "stats$accuracy is not a single numeric for ", result$gse,
      ": ", paste(utils::capture.output(str(result$stats)), collapse = " ")
    )
  }

  deviation <- if (!is.null(result$deviation) && !is.null(result$deviation$mean)) {
    as.numeric(result$deviation$mean)
  } else {
    NULL
  }

  rows <- list(
    list(
      step = "overall",
      cellsSynthesized = as.integer(result$n_synth_cells),
      accuracy = accuracy,
      deviation = deviation
    )
  )

  # Per-cluster synthesized counts
  cc <- result$cluster_counts
  if (length(cc) > 0) {
    cluster_names <- names(cc)
    if (is.null(cluster_names)) {
      cluster_names <- as.character(seq_along(cc))
    }
    for (i in seq_along(cc)) {
      rows[[length(rows) + 1]] <- list(
        step = paste0("cluster_", cluster_names[[i]]),
        cellsSynthesized = as.integer(cc[[i]]),
        accuracy = NULL,
        deviation = NULL
      )
    }
  }

  rows
}

export_one <- function(ds) {
  rds_path <- results_rds_path(ds)
  if (!file.exists(rds_path)) {
    stop("Results RDS not found: ", rds_path, "\nRun 03_run_scgft.R first.")
  }

  message("\nExporting JSON for ", ds$gse, " from ", rds_path)
  result <- readRDS(rds_path)
  rows <- build_export_rows(result)

  out <- results_json_path(ds)
  jsonlite::write_json(
    rows,
    path = out,
    pretty = TRUE,
    auto_unbox = TRUE,
    null = "null"
  )
  message("Wrote: ", out)
  message(jsonlite::toJSON(rows, pretty = TRUE, auto_unbox = TRUE, null = "null"))
  invisible(out)
}

main <- function() {
  message("Stage 04: export JSON")
  for (ds in DATASETS) {
    export_one(ds)
  }
  message("\nStage 04 complete.")
}

main()
