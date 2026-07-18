# scGFT Augmentation Pipeline

Standalone R pipeline that downloads two public scRNA-seq / snRNA-seq datasets from GEO, runs a standard Seurat preprocessing workflow, synthesizes additional cells with [scGFT](https://github.com/Sanofi-Public/PMCB-scGFT) (`Sanofi-Public/PMCB-scGFT`), evaluates synthesis quality, and exports metrics as JSON for a separate frontend.

This repo is self-contained — it is not part of another project.

## Datasets

| GEO accession | Disease / context | Assay |
|---|---|---|
| **GSE165577** | Rett syndrome (MECP2 mutant); cortex + ganglionic eminence organoids | scRNA-seq |
| **GSE156498** | Duchenne muscular dystrophy (DMD mutant vs wild-type); mouse TA muscle | snRNA-seq |

Supplementary files are pulled with `GEOquery::getGEOSuppFiles()`. Formats are **not assumed** up front: stage `01` downloads and prints a file listing; stage `02` detects whether each GSE is a 10x mtx triplet, an `.h5`, or a flat count matrix before loading.

## Pipeline stages

| Script | Role |
|---|---|
| `R/00_config.R` | Accessions, disease labels, paths (`data/`, `results/`), QC + scGFT parameters |
| `R/01_download_data.R` | Download GEO supplementary files; print / save listings |
| `R/02_preprocess.R` | Format-aware load → Seurat QC + Normalize / HVG / Scale / PCA / neighbors / clusters / UMAP → `data/processed/*.rds` |
| `R/03_run_scgft.R` | Install scGFT if needed → `RunScGFT()` → re-cluster → `statsScGFT()` → `results/*_scgft.rds` |
| `R/04_export_json.R` | Write per-dataset JSON arrays of `{step, cellsSynthesized, accuracy, deviation}` |

Metrics are taken from real package output only. `statsScGFT()` returns accuracy; deviation is parsed from `RunScGFT()` console messages (the package prints deviation but does not return it). Missing or unexpected structure is raised as an error — nothing is filled with placeholder numbers.

## Requirements

- R (≥ 4.x recommended) with `Rscript` on `PATH`
- Packages:
  - CRAN: `Seurat`, `SeuratObject`, `Matrix`, `jsonlite`, `devtools`
  - Bioconductor: `GEOquery`
  - GitHub: `Sanofi-Public/PMCB-scGFT` (installed automatically by stage `03` if missing)

Optional but useful: `data.table` (faster flat-matrix reads).

Example install:

```r
install.packages(c("Seurat", "SeuratObject", "Matrix", "jsonlite", "devtools", "data.table"))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("GEOquery")
```

## How to run

From the repo root (Linux / macOS / WSL):

```bash
chmod +x run_pipeline.sh
./run_pipeline.sh
```

The shell script runs the four stages in order with `Rscript`, stops on the first error, and prints which stage failed.

You can also run stages individually (useful after inspecting download listings):

```bash
Rscript R/01_download_data.R
# review data/raw/<GSE>/file_listing.txt, then:
Rscript R/02_preprocess.R
Rscript R/03_run_scgft.R
Rscript R/04_export_json.R
```

## Outputs

- `data/raw/<GSE>/` — GEO supplementary files + `file_listing.txt`
- `data/processed/<GSE>_<disease>.rds` — processed Seurat objects
- `results/<GSE>_<disease>_scgft.rds` — synthesis object + stats
- `results/<GSE>_<disease>.json` — frontend-facing metrics

Example JSON shape:

```json
[
  {
    "step": "overall",
    "cellsSynthesized": 1234,
    "accuracy": 94.03,
    "deviation": 0.15
  },
  {
    "step": "cluster_0",
    "cellsSynthesized": 120,
    "accuracy": null,
    "deviation": null
  }
]
```

Per-cluster rows include synthesized cell counts only; global `accuracy` / `deviation` live on the `overall` step because that is what scGFT reports.

## Configuration knobs

Edit `R/00_config.R` for QC thresholds (`QC_*`), PCA dims, cluster resolution, and synthesis scale (`SYNTH_SCALE`, `SYNTH_NCPMNTS`). Full 1× synthesis on large datasets can be slow and memory-heavy; lower `SYNTH_SCALE` for smoke tests.

## Notes

- `data/` and `results/` are gitignored.
- Stage `03` re-runs variable features / PCA / neighbors / clustering on the combined real+synthetic object before `statsScGFT()`, matching the scGFT README evaluation pattern.
