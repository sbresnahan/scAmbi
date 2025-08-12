# scAmbi — Mapping-ambiguity overdispersion correction for scRNA-seq

**scAmbi** provides tools to estimate and correct **overdispersion** arising from read-to-transcript mapping ambiguity (e.g., Salmon Alevin cell-level bootstraps), then evaluate effects on variability using within- and between-sample BCV analyses. It includes helpers to build corrected Seurat objects and plots for quick diagnostics.

---

## Why scAmbi?

Multi-mapping fragments and assignment ambiguity inflate inferential variance. Alevin can emit per-cell bootstrap replicates that capture this. scAmbi:

- computes per-gene OD from bootstraps (sparse-aware, block-wise),
- integrates bootstrap-, moments-, and prior-based estimates,
- constructs a corrected assay (counts scaled by 1/OD),
- quantifies improvement via edgeR BCV (within/between sample),
- offers plotting + summaries for rapid QC.

---

## Key features

- **Bootstrap OD (sparse-aware):** `compute_overdisp_sparse_aware()`
- **Integrated OD estimator:** `estimate_overdispersion_integrated()`
- **Seurat integration:** `process_and_create_seurat_corrected_improved()`
- **Within-sample BCV:** `calculate_within_sample_bcv()`, `analyze_within_sample_bcv()`
- **Between-sample BCV:** `extract_and_pseudobulk()`, `calculate_bcv_direct()`
- **Visualization:** `plot_within_sample_bcv()`, `plot_within_sample_summary()`, `plot_bcv()`, `plot_bcv_comparison()`
- **Utilities:** `read_eds_gc()`, `read_sample_data_improved()`, `set_feature_metadata()`, `extract_feature_vector()`

---

## Installation

### From GitHub (recommended)

```r
# install.packages("remotes")
remotes::install_github("sbresnahan/scAmbi")
```

### From a local source tarball/zip

```r
# If you have scAmbi.zip or a source tar.gz:
install.packages("scAmbi.zip", repos = NULL, type = "source")
# or use devtools:
# devtools::install("path/to/scAmbi/")
```

### Requirements

This package requires **R version 4.2 or higher**.

#### Core Dependencies

The following packages are automatically installed when you install `scAmbi`.

* `Seurat`
* `SeuratObject`
* `edgeR`
* `Matrix`
* `ggplot2`
* `patchwork`
* `eds`
* `dplyr`
* `tidyr`
* `magrittr`
* `parallel`
* `Rcpp`
* `rtracklayer`

#### Development Dependencies

These packages are necessary for building the vignettes and running tests. You can install them by running the following command:

```R
install.packages(c("knitr", "rmarkdown", "testthat"))

---

## Getting started

```r
library(scAmbi)

# 1) Estimate integrated overdispersion from one Alevin sample
#    (requires bootstraps in <sample>/alevin/quants_boot_mat.gz)
alevin_dir <- "path/to/<sample>/alevin"
# counts <- ... (genes x cells dgCMatrix), feats <- rownames(counts), cells <- colnames(counts)
od <- estimate_overdispersion_integrated(
  counts     = counts,
  alevin_dir = alevin_dir,
  n_boot     = 20,
  n_cores    = 4
)

# 2) Build a Seurat object with corrected assay (RNA_corr = counts / OD)
seu <- process_and_create_seurat_corrected_improved(
  sample_id = "S1", counts = counts, od = od, feats = feats, cells = cells
)

# 3) Within-sample BCV comparison (raw vs corrected)
wres <- analyze_within_sample_bcv(list(S1 = seu), assay_names = c("RNA", "RNA_corr"), n_groups = 10)
p <- plot_within_sample_bcv(wres, sample_name = "S1")
print(p)
```

[See the vignette](https://seantbresnahan.com/scambi) for a full, reproducible walkthrough.

---

## Alevin settings (inputs expected)

To use bootstrap-based OD, run Alevin with cell-level bootstraps enabled so that `alevin/quants_boot_mat.gz` is present:

```bash
salmon alevin \
  -l ISR \
  -1 <R1.fastq.gz> -2 <R2.fastq.gz> \
  --chromiumV3 \
  -i <txindex> \
  --whitelist <whitelist.txt> \
  --numCellBootstraps 20 \
  --dumpFeatures
```

Notes:
- scAmbi reads the **boot matrix** and associated index files via `eds::readEDS()`.
- For transcript-centric work, provide a suitable index/mapping to Alevin.

---

## Vignette

```r
# After install:
browseVignettes("scAmbi")
# Or build from source:
devtools::build_vignettes(); browseVignettes("scAmbi")
```

[The vignette](https://seantbresnahan.com/scambi) demonstrates OD estimation, Seurat correction, and BCV diagnostics end-to-end.

---

## License

GPL-3.  
**Maintainer:** Sean T. Bresnahan <stbresnahan@mdanderson.org>.

---

## Citation

If you use scAmbi, please cite this repository and the tools it builds upon (e.g., Salmon/Alevin, edgeR). A formal citation will be added once a preprint is available.
