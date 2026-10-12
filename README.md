# *Plasmodium falciparum* persister scRNA-seq — guided QC and exploration

This repository contains **one R script** for ten Cell Ranger gene-expression libraries (`6_1`, `6_2`, `7_1`, `7_2`, `8_1`, `8_2`, `9_1`, `9_2`, `10_1`, `10_2`). The analysis is **exploratory** and does not establish a persister-specific transcriptional state.

## Run in RStudio

Install dependencies once:

```r
install.packages(c('Seurat','SeuratObject','hdf5r','Matrix','ggplot2'))
```

Open `R/Pf_persister_guided_pipeline.R` in RStudio and click **Source**, or run:

```r
source('R/Pf_persister_guided_pipeline.R')
```

Change `data_root` in the script if needed (default `~/OneDrive/scRNA_data`). Each sample folder must contain `filtered_feature_bc_matrix.h5`. Results are written to `Pf_persister_guided_results` under the same root.

Plots are explicitly `print()`ed to the **RStudio Plots pane** (typically lower right) during an interactive run, and exported as PDF and PNG. Use the Plots pane's history arrows to inspect earlier figures. If your workspace pane layout differs, change it under **Tools → Global Options → Pane Layout**. `pause_after_plot <- TRUE` optionally prompts between figures, but interactive prompts may depend on how the script is sourced.

## Stages, purpose, expected results

| Stage | Why | Expected output |
|---|---|---|
| 1 | Verify installed dependencies, sample manifest | `sample_manifest_TO_ANNOTATE.csv` |
| 2 | Validate count matrices, barcodes and gene complexity | Raw metrics CSV and individual scatter plots |
| 3 | Identify possible mitochondrial annotations | Full feature audit, candidates, name mapping |
| 4 | Optionally remove verified mitochondrial genes | Before/after QC CSVs, Seurat RDS, QC figures |
| 5 | Avoid arbitrary 200-feature filtering of low-RNA cells | QC-ready Seurat objects and status log |
| 6 | Explore structure within each library | Per-library normalized RDS, PCA and UMAP |
| 7 | Compare raw QC across all libraries | Combined figures, diagnostics, pipeline log |

**Mitochondrial annotation:** `^MT-` previously matched zero genes. Candidate pattern matches are exported for review, not automatically excluded. Set `verified_mt_file` only after confirming gene identities using your actual reference. The file may list original HDF5 feature names or their Seurat equivalents. If no verified list is provided, mitochondrial genes are **not removed**, and `percent.mt` is not invented.

**Cell filtering:** no threshold is applied by default. To filter, specify `filter_thresholds_file` with columns `sample_id,min_features,max_features,min_umis,max_umis,max_percent_mt`; blanks indicate no cutoff. Sample `8_2` previously had low gene detection and merits specific review before filtering.

**Experimental metadata:** complete the manifest using actual treatment groups, replicate relationships, timepoints and batches. Folder names alone cannot establish comparisons. No cross-sample integration, cluster identification, or differential expression is included.

**Privacy:** Do not commit patient-level or unpublished sequencing data without authorization. `*.h5`, `*.rds` and results directories are excluded by `.gitignore`.

**Validation:** Code has been reviewed for consistency but **not executed** against the ten local HDF5 matrices. Run and review `pipeline_status.csv` and QC plots before downstream work.
