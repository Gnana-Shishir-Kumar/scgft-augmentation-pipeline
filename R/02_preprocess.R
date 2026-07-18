# Load GEO supplementary data using format detection, run Seurat preprocessing,
# and save processed objects to data/processed/.

source(file.path(
  if (basename(getwd()) == "R") ".." else ".",
  "R", "00_config.R"
))

required_pkgs <- c("Seurat", "SeuratObject", "Matrix")
missing <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) {
  stop(
    "Missing required packages: ", paste(missing, collapse = ", "), "\n",
    "Install with install.packages(c('Seurat', 'SeuratObject', 'Matrix'))."
  )
}

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
})

# ---------------------------------------------------------------------------
# Format detection (based on what 01_download_data.R actually retrieved)
# ---------------------------------------------------------------------------

is_10x_dir <- function(dir_path) {
  files <- list.files(dir_path, full.names = FALSE)
  files_l <- tolower(files)
  has_mtx <- any(grepl("matrix\\.mtx(\\.gz)?$", files_l))
  has_barcodes <- any(grepl("barcodes\\.tsv(\\.gz)?$", files_l))
  has_features <- any(grepl("(features|genes)\\.tsv(\\.gz)?$", files_l))
  has_mtx && has_barcodes && has_features
}

find_10x_dirs <- function(root_dir) {
  # Include root and all subdirectories
  dirs <- unique(c(
    root_dir,
    list.dirs(root_dir, recursive = TRUE, full.names = TRUE)
  ))
  dirs[vapply(dirs, is_10x_dir, logical(1))]
}

find_h5_files <- function(root_dir) {
  list.files(
    root_dir,
    pattern = "\\.(h5|hdf5)$",
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = TRUE
  )
}

find_flat_matrices <- function(root_dir) {
  # Prefer uncompressed / gz text matrices; skip tiny readmes
  candidates <- list.files(
    root_dir,
    pattern = "\\.(txt|tsv|csv|txt\\.gz|tsv\\.gz|csv\\.gz)$",
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = TRUE
  )
  # Drop GEO soft/readme-ish names when possible
  skip <- grepl("(readme|file_listing|soft|series_matrix)", basename(candidates), ignore.case = TRUE)
  candidates <- candidates[!skip]
  # Prefer names that look like counts / expression
  prefer <- grepl("(count|express|matrix|umi|gene)", basename(candidates), ignore.case = TRUE)
  if (any(prefer)) {
    candidates <- c(candidates[prefer], candidates[!prefer])
    candidates <- unique(candidates)
  }
  candidates
}

detect_format <- function(root_dir) {
  tenx <- find_10x_dirs(root_dir)
  h5   <- find_h5_files(root_dir)
  flat <- find_flat_matrices(root_dir)

  message("Format detection under: ", root_dir)
  message("  10x directories: ", if (length(tenx)) paste(tenx, collapse = " | ") else "(none)")
  message("  H5 files:     ", if (length(h5)) paste(h5, collapse = " | ") else "(none)")
  message("  Flat matrices:", if (length(flat)) paste(flat, collapse = " | ") else "(none)")

  if (length(tenx) > 0) {
    list(type = "10x", paths = tenx)
  } else if (length(h5) > 0) {
    list(type = "h5", paths = h5)
  } else if (length(flat) > 0) {
    list(type = "matrix", paths = flat)
  } else {
    list(type = "unknown", paths = character())
  }
}

# ---------------------------------------------------------------------------
# Loaders
# ---------------------------------------------------------------------------

load_10x <- function(dirs) {
  objs <- lapply(seq_along(dirs), function(i) {
    d <- dirs[[i]]
    message("Reading 10X data from: ", d)
    counts <- Seurat::Read10X(data.dir = d)
    if (is.list(counts) && !is.data.frame(counts)) {
      # Multi-modal 10x; prefer Gene Expression
      if ("Gene Expression" %in% names(counts)) {
        counts <- counts[["Gene Expression"]]
      } else {
        counts <- counts[[1]]
      }
    }
    project <- paste0("s", i)
    Seurat::CreateSeuratObject(
      counts = counts,
      project = project,
      min.cells = SEURAT_MIN_CELLS,
      min.features = SEURAT_MIN_FEATURES
    )
  })
  if (length(objs) == 1) {
    return(objs[[1]])
  }
  message("Merging ", length(objs), " 10X samples...")
  merged <- merge(objs[[1]], y = objs[-1], add.cell.ids = paste0("s", seq_along(objs)))
  # Seurat v5 split layers need joining before a single-assay NormalizeData
  if ("JoinLayers" %in% getNamespaceExports("SeuratObject")) {
    merged <- SeuratObject::JoinLayers(merged)
  }
  merged
}

load_h5 <- function(files) {
  objs <- lapply(seq_along(files), function(i) {
    f <- files[[i]]
    message("Reading 10X H5 from: ", f)
    counts <- Seurat::Read10X_h5(filename = f)
    if (is.list(counts) && !is.data.frame(counts)) {
      if ("Gene Expression" %in% names(counts)) {
        counts <- counts[["Gene Expression"]]
      } else {
        counts <- counts[[1]]
      }
    }
    Seurat::CreateSeuratObject(
      counts = counts,
      project = paste0("h5", i),
      min.cells = SEURAT_MIN_CELLS,
      min.features = SEURAT_MIN_FEATURES
    )
  })
  if (length(objs) == 1) {
    return(objs[[1]])
  }
  message("Merging ", length(objs), " H5 samples...")
  merged <- merge(objs[[1]], y = objs[-1], add.cell.ids = paste0("h5", seq_along(objs)))
  if ("JoinLayers" %in% getNamespaceExports("SeuratObject")) {
    merged <- SeuratObject::JoinLayers(merged)
  }
  merged
}

read_table_auto <- function(path) {
  message("Reading flat matrix: ", path)
  sep <- if (grepl("\\.csv(\\.gz)?$", path, ignore.case = TRUE)) "," else "\t"
  # data.table is faster when available; base fallback otherwise
  if (requireNamespace("data.table", quietly = TRUE)) {
    dt <- data.table::fread(path, sep = sep, header = TRUE, data.table = FALSE)
  } else {
    gz <- grepl("\\.gz$", path, ignore.case = TRUE)
    con <- if (gz) gzfile(path, open = "rt") else file(path, open = "rt")
    on.exit(close(con), add = TRUE)
    dt <- utils::read.table(con, sep = sep, header = TRUE, check.names = FALSE, quote = "")
  }
  # First column as gene names when non-numeric
  if (!is.numeric(dt[[1]])) {
    rn <- as.character(dt[[1]])
    dt <- dt[, -1, drop = FALSE]
    mat <- as.matrix(dt)
    storage.mode(mat) <- "numeric"
    rownames(mat) <- make.unique(rn)
  } else {
    mat <- as.matrix(dt)
    storage.mode(mat) <- "numeric"
  }
  mat
}

load_flat_matrix <- function(files) {
  # Use the first plausible matrix; warn if several
  if (length(files) > 1) {
    message(
      "Multiple flat matrices found; using the first preferred candidate:\n  ",
      files[[1]],
      "\nOther candidates:\n  ",
      paste(files[-1], collapse = "\n  ")
    )
  }
  mat <- read_table_auto(files[[1]])

  # Heuristic: genes x cells is typical; if more rows than cols looks like cells x genes, transpose
  if (nrow(mat) < ncol(mat) && nrow(mat) < 500) {
    message(
      "Matrix has ", nrow(mat), " rows x ", ncol(mat),
      " cols — fewer rows than expected for genes; transposing to genes x cells."
    )
    mat <- t(mat)
  }

  Seurat::CreateSeuratObject(
    counts = Matrix::Matrix(mat, sparse = TRUE),
    project = "matrix",
    min.cells = SEURAT_MIN_CELLS,
    min.features = SEURAT_MIN_FEATURES
  )
}

load_dataset <- function(gse) {
  root <- raw_gse_dir(gse)
  if (!dir.exists(root)) {
    stop("Raw data directory missing for ", gse, ": ", root, "\nRun 01_download_data.R first.")
  }

  fmt <- detect_format(root)
  if (identical(fmt$type, "unknown")) {
    stop(
      "Could not detect a loadable expression format under ", root, ".\n",
      "Re-run 01_download_data.R and inspect file_listing.txt."
    )
  }

  message("Using detected format: ", fmt$type)
  switch(
    fmt$type,
    "10x"    = load_10x(fmt$paths),
    "h5"     = load_h5(fmt$paths),
    "matrix" = load_flat_matrix(fmt$paths),
    stop("Unhandled format: ", fmt$type)
  )
}

# ---------------------------------------------------------------------------
# Seurat pipeline + QC
# ---------------------------------------------------------------------------

run_qc_and_preprocess <- function(sobj, ds) {
  message("\n--- Preprocess ", ds$gse, " (", ds$disease, ") ---")
  message("Cells before QC: ", ncol(sobj), " | Features: ", nrow(sobj))

  sobj[["percent.mt"]] <- Seurat::PercentageFeatureSet(sobj, pattern = "^MT-|^mt-|^Mt-")

  message(
    "QC filters: nFeature_RNA in [", QC_MIN_FEATURES, ", ", QC_MAX_FEATURES,
    "], percent.mt < ", QC_MAX_MT
  )
  sobj <- subset(
    sobj,
    subset = nFeature_RNA >= QC_MIN_FEATURES &
      nFeature_RNA <= QC_MAX_FEATURES &
      percent.mt < QC_MAX_MT
  )
  message("Cells after QC: ", ncol(sobj), " | Features: ", nrow(sobj))

  if (ncol(sobj) < 50) {
    stop(
      "Too few cells after QC for ", ds$gse, " (n=", ncol(sobj), "). ",
      "Adjust QC_* in 00_config.R or inspect the loaded matrix."
    )
  }

  sobj <- Seurat::NormalizeData(sobj, normalization.method = "LogNormalize", scale.factor = 1e4)
  sobj <- Seurat::FindVariableFeatures(sobj, selection.method = "vst", nfeatures = 2000)
  sobj <- Seurat::ScaleData(sobj, verbose = FALSE)
  sobj <- Seurat::RunPCA(sobj, npcs = max(PCA_DIMS), seed.use = UMAP_SEED, verbose = FALSE)
  sobj <- Seurat::FindNeighbors(sobj, dims = PCA_DIMS, verbose = FALSE)
  sobj <- Seurat::FindClusters(sobj, resolution = CLUSTER_RES, random.seed = UMAP_SEED, verbose = FALSE)
  sobj <- Seurat::RunUMAP(sobj, dims = PCA_DIMS, seed.use = UMAP_SEED, verbose = FALSE)

  message(
    "Clusters (seurat_clusters): ",
    paste(names(table(sobj$seurat_clusters)), table(sobj$seurat_clusters), sep = "=", collapse = ", ")
  )
  sobj
}

main <- function() {
  message("Stage 02: preprocess Seurat objects")

  for (ds in DATASETS) {
    message("\n########################################")
    message("Dataset: ", ds$gse, " — ", ds$note)
    message("########################################")

    sobj <- load_dataset(ds$gse)
    sobj <- run_qc_and_preprocess(sobj, ds)

    out <- processed_rds_path(ds)
    saveRDS(sobj, file = out)
    message("Saved processed Seurat object: ", out)
  }

  message("\nStage 02 complete.")
}

main()
