# Run scGFT synthesis + quality stats on each processed Seurat object.
# Does not invent metrics: accuracy comes from statsScGFT(); deviation is
# parsed from RunScGFT() messages (the package prints it but does not return it).

source(file.path(
  if (basename(getwd()) == "R") ".." else ".",
  "R", "00_config.R"
))

ensure_scgft <- function() {
  if (!requireNamespace("scgft", quietly = TRUE)) {
    message("Package 'scgft' not installed — installing from Sanofi-Public/PMCB-scGFT ...")
    if (!requireNamespace("devtools", quietly = TRUE)) {
      install.packages("devtools", repos = "https://cloud.r-project.org")
    }
    devtools::install_github(
      "Sanofi-Public/PMCB-scGFT",
      build_vignettes = FALSE,
      upgrade = "never"
    )
  }
  if (!requireNamespace("scgft", quietly = TRUE)) {
    stop("Failed to install/load scgft from Sanofi-Public/PMCB-scGFT")
  }
  suppressPackageStartupMessages(library(scgft))
}

required_pkgs <- c("Seurat", "SeuratObject")
missing <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) {
  stop("Missing required packages: ", paste(missing, collapse = ", "))
}

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
})

ensure_scgft()

#' Capture "Deviation from originals: X +/- Y" messages emitted by RunScGFT.
parse_deviation_messages <- function(msgs) {
  hits <- grep("Deviation from originals:", msgs, value = TRUE, fixed = TRUE)
  if (length(hits) == 0) {
    return(NULL)
  }
  # Example: "Deviation from originals: 0.15 +/- 0.01"
  vals <- lapply(hits, function(m) {
    m <- sub(".*Deviation from originals:\\s*", "", m)
    parts <- strsplit(m, "\\+/-")[[1]]
    mean_v <- suppressWarnings(as.numeric(trimws(parts[1])))
    sd_v <- if (length(parts) >= 2) {
      suppressWarnings(as.numeric(trimws(parts[2])))
    } else {
      NA_real_
    }
    c(mean = mean_v, sd = sd_v)
  })
  means <- vapply(vals, `[[`, numeric(1), "mean")
  sds   <- vapply(vals, `[[`, numeric(1), "sd")
  if (anyNA(means)) {
    return(NULL)
  }
  list(
    mean = mean(means),
    sd = if (all(is.na(sds))) NA_real_ else mean(sds, na.rm = TRUE),
    raw_messages = hits
  )
}

run_scgft_one <- function(ds) {
  in_path <- processed_rds_path(ds)
  if (!file.exists(in_path)) {
    stop("Processed object not found: ", in_path, "\nRun 02_preprocess.R first.")
  }

  message("\n########################################")
  message("scGFT: ", ds$gse, " (", ds$disease, ")")
  message("Loading: ", in_path)
  message("########################################")

  sobj <- readRDS(in_path)
  n_cells <- ncol(sobj)
  nsynth <- as.integer(round(SYNTH_SCALE * n_cells))
  if (nsynth < 1) {
    stop("nsynth computed as ", nsynth, " — increase SYNTH_SCALE in 00_config.R")
  }

  message(
    "Running RunScGFT(nsynth=", format(nsynth, big.mark = ","),
    ", ncpmnts=", SYNTH_NCPMNTS,
    ", groups='", SYNTH_GROUPS, "') on ",
    format(n_cells, big.mark = ","), " cells ..."
  )
  message("(scGFT prints synthesis progress / deviation as it runs)")

  if (!SYNTH_GROUPS %in% colnames(sobj[[]])) {
    stop("Metadata column '", SYNTH_GROUPS, "' missing from processed object.")
  }

  # Collect messages so we can recover deviation (printed, not returned)
  msg_buf <- character()
  sobj_synt <- withCallingHandlers(
    {
      scgft::RunScGFT(
        object  = sobj,
        nsynth  = nsynth,
        ncpmnts = SYNTH_NCPMNTS,
        groups  = SYNTH_GROUPS
      )
    },
    message = function(m) {
      msg_buf <<- c(msg_buf, conditionMessage(m))
    }
  )

  deviation <- parse_deviation_messages(msg_buf)
  if (is.null(deviation)) {
    warning(
      "Could not parse deviation from RunScGFT messages for ", ds$gse, ". ",
      "Deviation will be stored as NULL (not fabricated)."
    )
  } else {
    message(
      "Parsed deviation mean=", round(deviation$mean, 4),
      " sd=", round(deviation$sd, 4)
    )
  }

  n_synth_cells <- sum(sobj_synt$synthesized == "yes", na.rm = TRUE)
  message("Synthesized cell count in object: ", format(n_synth_cells, big.mark = ","))

  # Re-embed / re-cluster combined real+synth object (scGFT README pattern)
  # so statsScGFT compares post-synthesis cluster co-membership.
  message("Re-running variable features / PCA / neighbors / clusters on combined object...")
  sobj_synt <- Seurat::FindVariableFeatures(sobj_synt, selection.method = "vst", nfeatures = 2000)
  sobj_synt <- Seurat::ScaleData(sobj_synt, verbose = FALSE)
  sobj_synt <- Seurat::RunPCA(sobj_synt, npcs = max(PCA_DIMS), seed.use = UMAP_SEED, verbose = FALSE)
  sobj_synt <- Seurat::FindNeighbors(sobj_synt, dims = PCA_DIMS, verbose = FALSE)
  sobj_synt <- Seurat::FindClusters(
    sobj_synt,
    resolution = CLUSTER_RES,
    random.seed = UMAP_SEED,
    verbose = FALSE
  )
  sobj_synt <- Seurat::RunUMAP(sobj_synt, dims = PCA_DIMS, seed.use = UMAP_SEED, verbose = FALSE)

  message("Running statsScGFT(groups='", SYNTH_GROUPS, "') ...")
  stats <- tryCatch(
    scgft::statsScGFT(object = sobj_synt, groups = SYNTH_GROUPS),
    error = function(e) {
      stop(
        "statsScGFT() failed for ", ds$gse, ": ", conditionMessage(e),
        call. = FALSE
      )
    }
  )

  # Surface unexpected return shapes instead of inventing fields
  if (!is.list(stats) || is.null(stats$accuracy)) {
    stop(
      "Unexpected statsScGFT() return for ", ds$gse, ".\n",
      "Expected list with $accuracy; got: ",
      paste(utils::capture.output(str(stats)), collapse = "\n")
    )
  }
  message("statsScGFT accuracy (%): ", stats$accuracy)

  # Per-cluster synthesized cell counts (from metadata after synthesis)
  synth_only <- subset(sobj_synt, subset = synthesized == "yes")
  cluster_counts <- as.list(table(synth_only[[SYNTH_GROUPS]]))

  result <- list(
    gse = ds$gse,
    disease = ds$disease,
    label = ds$label,
    n_original_cells = n_cells,
    n_synth_requested = nsynth,
    n_synth_cells = n_synth_cells,
    ncpmnts = SYNTH_NCPMNTS,
    groups = SYNTH_GROUPS,
    deviation = deviation,          # list(mean, sd, raw_messages) or NULL
    stats = stats,                  # list(accuracy=...) from package
    cluster_counts = cluster_counts,
    sobj_synt = sobj_synt
  )

  out <- results_rds_path(ds)
  saveRDS(result, file = out)
  message("Saved scGFT results: ", out)
  invisible(result)
}

main <- function() {
  message("Stage 03: RunScGFT + statsScGFT")
  for (ds in DATASETS) {
    run_scgft_one(ds)
  }
  message("\nStage 03 complete.")
}

main()
