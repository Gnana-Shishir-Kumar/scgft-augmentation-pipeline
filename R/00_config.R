# Project configuration for the scGFT augmentation pipeline.
# Sourced by the other R/*.R scripts; not meant to be run alone.

options(stringsAsFactors = FALSE)
set.seed(42)

# Resolve project root whether launched from repo root or from R/
PROJECT_ROOT <- local({
  wd <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  if (basename(wd) == "R") {
    normalizePath(dirname(wd), winslash = "/", mustWork = TRUE)
  } else {
    wd
  }
})

DATA_DIR      <- file.path(PROJECT_ROOT, "data")
RAW_DIR       <- file.path(DATA_DIR, "raw")
PROCESSED_DIR <- file.path(DATA_DIR, "processed")
RESULTS_DIR   <- file.path(PROJECT_ROOT, "results")

for (d in c(DATA_DIR, RAW_DIR, PROCESSED_DIR, RESULTS_DIR)) {
  if (!dir.exists(d)) {
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
    message("Created directory: ", d)
  }
}

# GEO accessions and disease labels for downstream JSON naming
DATASETS <- list(
  list(
    gse     = "GSE165577",
    disease = "Rett_syndrome",
    label   = "rett",
    note    = "MECP2 mutant; cortex + ganglionic eminence organoids (scRNA-seq)"
  ),
  list(
    gse     = "GSE156498",
    disease = "Duchenne_muscular_dystrophy",
    label   = "dmd",
    note    = "DMD mutant vs wild-type; mouse TA muscle (snRNA-seq)"
  )
)

# Seurat object creation / QC defaults (printed at runtime)
SEURAT_MIN_CELLS    <- 3L
SEURAT_MIN_FEATURES <- 200L

QC_MIN_FEATURES <- 200L
QC_MAX_FEATURES <- 6000L
QC_MAX_MT       <- 20  # percent mitochondrial counts

# Dimensionality reduction / clustering
PCA_DIMS       <- 1:30
CLUSTER_RES    <- 0.5
UMAP_SEED      <- 42L

# scGFT synthesis controls
# nsynth = SYNTH_SCALE * n_cells (1 = synthesize one synthetic cell per real cell)
SYNTH_SCALE <- 1
# Number of Fourier complex components to modify (see scGFT README examples)
SYNTH_NCPMNTS <- 10L
SYNTH_GROUPS  <- "seurat_clusters"

#' Build a stable basename for a dataset: GSE165577_Rett_syndrome
dataset_basename <- function(ds) {
  paste(ds$gse, ds$disease, sep = "_")
}

#' Path helpers
raw_gse_dir <- function(gse) file.path(RAW_DIR, gse)
processed_rds_path <- function(ds) {
  file.path(PROCESSED_DIR, paste0(dataset_basename(ds), ".rds"))
}
results_rds_path <- function(ds) {
  file.path(RESULTS_DIR, paste0(dataset_basename(ds), "_scgft.rds"))
}
results_json_path <- function(ds) {
  file.path(RESULTS_DIR, paste0(dataset_basename(ds), ".json"))
}

message("Project root: ", PROJECT_ROOT)
message("Datasets: ", paste(vapply(DATASETS, function(d) d$gse, character(1)), collapse = ", "))
