# One-shot installer for pipeline dependencies (CRAN + Bioconductor + scGFT).

repos <- "https://cloud.r-project.org"

# Prefer a user-writable library (system site-library often needs sudo)
user_lib <- Sys.getenv("R_LIBS_USER")
if (!nzchar(user_lib)) {
  user_lib <- path.expand(file.path("~", "R", paste0(R.version$platform, "-library"), paste(R.version$major, strsplit(R.version$minor, "[.]")[[1]][1], sep = ".")))
}
if (!dir.exists(user_lib)) {
  dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)
  message("Created user library: ", user_lib)
}
.libPaths(c(user_lib, .libPaths()))
message("Installing into: ", .libPaths()[[1]])

cran_pkgs <- c(
  "devtools",
  "Seurat",
  "SeuratObject",
  "Matrix",
  "jsonlite",
  "data.table"
)

message("Installing CRAN packages...")
for (pkg in cran_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message("Installing ", pkg, " ...")
    install.packages(pkg, repos = repos)
  } else {
    message(pkg, " already installed")
  }
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  message("Installing BiocManager ...")
  install.packages("BiocManager", repos = repos)
}

message("Installing GEOquery from Bioconductor ...")
BiocManager::install("GEOquery", ask = FALSE, update = FALSE)

if (!requireNamespace("scgft", quietly = TRUE)) {
  message("Installing scgft from GitHub (Sanofi-Public/PMCB-scGFT) ...")
  devtools::install_github(
    "Sanofi-Public/PMCB-scGFT",
    build_vignettes = FALSE,
    upgrade = "never"
  )
} else {
  message("scgft already installed")
}

needed <- c(cran_pkgs, "BiocManager", "GEOquery", "scgft")
ok <- vapply(needed, requireNamespace, logical(1), quietly = TRUE)
message("\nInstall check:")
print(ok)
if (!all(ok)) {
  stop("Missing packages: ", paste(needed[!ok], collapse = ", "))
}
message("\nAll pipeline dependencies are available.")
