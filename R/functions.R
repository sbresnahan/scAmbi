utils::globalVariables(c(
  "ave_log_cpm","bcv","bcv_trend","BCV","Median_BCV","Type","Assay",
  "RNA","RNA_corr","select","sample"
))

#' scAmbi: Overdispersion Correction for scRNA-seq Data
#'
#' @description
#' The `scAmbi` package provides a comprehensive toolkit for addressing technical
#' noise in single-cell RNA sequencing (scRNA-seq) data, specifically the
#' overdispersion caused by read-to-transcript mapping ambiguity. This issue is
#' common in workflows using tools like Salmon/Alevin, which employ bootstrap
#' replicates to quantify uncertainty.
#'
#' @details
#' Key functionalities include:
#' \itemize{
#'   \item \strong{Estimation}: Implements multiple methods to estimate per-gene overdispersion,
#'     including a fast, sparse-aware method for Alevin's bootstrap matrices, a
#'     method-of-moments approach, and a complexity-based prior derived from
#'     transcript annotations (GTF files).
#'   \item \strong{Integration}: Integrates these disparate estimates into a single, robust
#'     overdispersion factor for each transcript using a weighted geometric mean and
#'     empirical Bayes shrinkage.
#'   \item \strong{Correction}: Generates a new, corrected data assay within a Seurat object
#'     where raw counts are adjusted by the bootstrap scaling factors by default.
#'     Vectors for each estimate are also stored for manual correction.
#'   \item \strong{Evaluation}: Offers tools to calculate and visualize the Biological
#'     Coefficient of Variation (BCV) both within a single sample (across pseudo-bulk
#'     replicates) and between different samples. This allows for direct assessment of
#'     how the correction reduces technical variability.
#' }
#' The package is designed to integrate smoothly with Seurat-based analysis pipelines.
#'
#' @seealso For a detailed tutorial, see the package vignette:
#'   \code{vignette("scAmbi-intro", package = "scAmbi")}
#'
#' @docType package
#' @name scAmbi
#' @keywords internal
#'
#' @importFrom methods as
#' @importFrom stats cor median quantile sd var model.matrix
#' @importFrom utils read.table write.csv
#' @importFrom Matrix rowSums colSums rowMeans Diagonal
#' @importFrom Seurat CreateSeuratObject CreateAssayObject GetAssayData DefaultAssay DefaultAssay<-
#' @importFrom edgeR DGEList calcNormFactors filterByExpr estimateDisp
#' @importFrom ggplot2 ggplot aes geom_point geom_line geom_hline geom_abline geom_bar geom_boxplot geom_violin geom_text scale_color_gradient2 scale_fill_manual coord_cartesian coord_equal labs theme_classic theme_void theme element_text annotate ggsave
#' @import patchwork
#' @importFrom magrittr %>%
#' @importFrom tidyr pivot_wider
#' @importFrom dplyr select
#' @importFrom parallel mclapply detectCores
#' @importFrom rtracklayer import
#' @useDynLib scAmbi, .registration = TRUE
#' @importFrom Rcpp sourceCpp
"_PACKAGE"


#' Read an EDS matrix with expected dimensions and logging
#'
#' Reads a sparse EDS matrix from disk, checks that its dimensions match the
#' expected number of features and cells, and logs progress.
#'
#' @param path Character string; path to the `.eds` file to read.
#' @param n_feat Integer; expected number of features (rows) in the matrix.
#' @param n_cells Integer; expected number of cells (columns) in the matrix.
#' @param label Character string; label used in log messages to identify the dataset.
#'
#' @return A sparse matrix with dimensions \code{n_feat} x \code{n_cells}.
#' @importFrom eds readEDS
#' @examples
#' \dontrun{
#' m <- read_eds_gc("counts.eds", n_feat = 20000, n_cells = 5000, label = "Sample1")
#' }
#' @export


read_eds_gc <- function(path, n_feat, n_cells, label) {
  .log("[read] %s expecting features=%d cells=%d (file-native order)", label, n_feat, n_cells)
  m <- eds::readEDS(path, numOfGenes = n_feat, numOfOriginalCells = n_cells)  # features x cells
  .log("[read] %s got %d x %d (features x cells)", label, nrow(m), ncol(m))
  if (nrow(m) != n_feat || ncol(m) != n_cells) {
    stop(sprintf("[read] %s dimension mismatch: expected %dx%d, got %dx%d",
                 label, n_feat, n_cells, nrow(m), ncol(m)))
  }
  m
}


#' Calculate Gene Complexity from Transcript Annotations
#'
#' Computes the number of transcripts associated with each gene from a GTF or
#' a pre-processed transcript information file. This metric can be used to
#' inform a prior on mapping ambiguity.
#'
#' @param gtf_file Character string; path to a GTF file containing transcript
#'   annotations. This is used if `transcript_info_file` is not provided.
#' @param transcript_info_file Character string; path to a tab-delimited file
#'   with a `gene_id` column. If provided, this file is used in preference to
#'   the GTF file.
#'
#' @return A named numeric vector where each element is the number of transcripts
#'   for a gene and the name is the gene ID. Returns `NULL` if neither input
#'   file is provided or accessible.
#'
#' @importFrom rtracklayer import
#' @examples
#' \dontrun{
#' # From a pre-processed transcript info file
#' complexity <- get_gene_complexity(transcript_info_file = "tx_info.tsv")
#'
#' # From a GTF annotation file
#' complexity <- get_gene_complexity(gtf_file = "annotation.gtf")
#' }
#' @export


get_gene_complexity <- function(gtf_file = NULL, transcript_info_file = NULL) {
  if (!is.null(transcript_info_file) && file.exists(transcript_info_file)) {
    tx_info <- read.table(transcript_info_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    gene_complexity <- table(tx_info$gene_id)
    return(as.numeric(gene_complexity))
  } else if (!is.null(gtf_file) && file.exists(gtf_file)) {
    gtf <- rtracklayer::import(gtf_file)
    transcripts <- gtf[gtf$type == "transcript"]
    gene_complexity <- table(transcripts$gene_id)
    return(as.numeric(gene_complexity))
  } else {
    return(NULL)
  }
}



#' Compute Overdispersion from Alevin Bootstraps
#'
#' Estimates per-gene overdispersion from Alevin's bootstrap matrices using a
#' memory-efficient, block-wise computation strategy. This function averages
#' per-cell overdispersion estimates for each gene and applies light shrinkage
#' for robustness.
#'
#' @param alevin_dir Path to the Alevin output directory containing
#'   `quants_boot_mat.gz` and associated index files.
#' @param n_boot Integer; number of bootstrap replicates per cell.
#' @param block_cells Integer; number of cells to process in each parallel block
#'   (e.g., 256). A larger size may be faster but uses more memory.
#' @param n_cores Integer; number of parallel cores to use for processing.
#' @param omp_threads Integer; number of threads for low-level OpenMP/BLAS
#'   operations within the C++ code.
#' @param min_cells_expr Integer; minimum number of cells a gene must be
#'   expressed in to receive a non-trivial overdispersion estimate.
#' @param pseudocount Numeric; a small value added to the mean in the
#'   overdispersion calculation (`var / (mean + pseudocount)`) to prevent
#'   division by zero.
#' @param verbose Logical; if `TRUE`, print progress messages.
#' @param debug_first_block Logical; if `TRUE`, the first block is run in
#'   serial to provide clearer error messages for debugging.
#'
#' @return A list containing:
#' \describe{
#'   \item{`OverDisp`}{A numeric vector of the final per-gene overdispersion estimates.}
#'   \item{`DF`}{An integer vector of the number of cells contributing to each gene's estimate (degrees of freedom).}
#'   \item{`features`}{A character vector of the gene IDs.}
#'   \item{`expr_cells`}{A numeric vector of the total number of cells in which each gene was expressed.}
#' }
#' @export


compute_overdisp_sparse_aware <- function(
    alevin_dir,
    n_boot,
    block_cells    = 256L,
    n_cores        = max(1L, parallel::detectCores() - 1L),
    omp_threads    = 2L,
    min_cells_expr = 10,
    pseudocount    = 0.1,
    verbose        = TRUE,
    debug_first_block = FALSE
) {
  
  # --- files & dimensions ---
  f_boot <- file.path(alevin_dir, "quants_boot_mat.gz")
  f_cols <- file.path(alevin_dir, "quants_mat_cols.txt")
  f_rows <- file.path(alevin_dir, "quants_mat_rows.txt")
  stopifnot(file.exists(f_boot), file.exists(f_cols), file.exists(f_rows))
  
  genes <- scan(f_cols, what = "character", quiet = TRUE)
  cells <- scan(f_rows, what = "character", quiet = TRUE)
  G <- length(genes); C <- length(cells)
  .log("[OD] Expecting boot with C=%d, G=%d, N=%d", C, G, n_boot)
  
  # Read the full bootstrap matrix
  .log("[OD] Reading the full boot matrix into memory.")
  boot <- eds::readEDS(f_boot, numOfGenes = G, numOfOriginalCells = C * n_boot)
  if (!inherits(boot, "dgCMatrix")) boot <- as(boot, "dgCMatrix")
  .log("[OD] boot dims = %d x %d (genes x C*N)", nrow(boot), ncol(boot))
  
  # Sanity on shape
  if (ncol(boot) != C * n_boot)
    stop(sprintf("[OD] ncol(boot)=%d != C*n_boot=%d", ncol(boot), C * n_boot))
  
  # Setup for block processing
  Sys.setenv(OMP_NUM_THREADS = as.integer(omp_threads))
  
  # Create blocks of 0-based cell indices [start_cell0, end_cell0]
  starts <- seq.int(1L, C, by = block_cells)
  blocks <- vector("list", length(starts))
  for (i in seq_along(starts)) {
    s <- starts[i]
    e <- min(C, s + block_cells - 1L)
    blocks[[i]] <- c(s - 1L, e - 1L)
  }
  
  block_fn <- function(bl) {
    start_col <- bl[1] * n_boot
    end_col   <- (bl[2] + 1L) * n_boot   # exclusive
    compute_overdisp_block_improved(boot, start_col, end_col, n_boot, G, 1e-8, pseudocount)
  }
  
  # Optional: run the first block single-core for clear error messages
  if (debug_first_block) {
    bl <- blocks[[1]]
    start_col <- bl[1] * n_boot
    end_col   <- (bl[2] + 1L) * n_boot
    .log("[debug] first block: start_cell0=%d end_cell0=%d | start_col=%d end_col(excl)=%d of [0,%d)",
         bl[1], bl[2], start_col, end_col, ncol(boot))
    invisible(block_fn(bl))  # will stop on error with a clear message
  }
  
  block_fn_safe <- function(bl) {
    tryCatch(
      block_fn(bl),
      error = function(e) {
        stop(sprintf("[block cells %d..%d] %s", bl[1], bl[2], conditionMessage(e)), call. = FALSE)
      }
    )
  }
  
  # Parallel processing
  if (n_cores > 1) {
    out_list <- parallel::mclapply(
      blocks, block_fn_safe,
      mc.cores = n_cores,
      mc.preschedule = FALSE,
      mc.silent = FALSE,
      mc.set.seed = FALSE
    )
  } else {
    out_list <- lapply(blocks, block_fn_safe)
  }
  
  # Combine results
  OverDisp_sum <- Reduce(`+`, lapply(out_list, `[[`, "OD_sum"))
  DF          <- Reduce(`+`, lapply(out_list, `[[`, "DF_sum"))
  expr_cells  <- Reduce(`+`, lapply(out_list, `[[`, "expr_cells"))
  
  .log("[OD] Processed all %d cells", C)
  
  # Compute overdispersion with expression filtering
  OverDisp <- rep_len(1.0, G)
  expressed <- which(expr_cells >= min_cells_expr)
  
  if (length(expressed) > 0) {
    od_expressed <- OverDisp_sum[expressed] / DF[expressed]
    OverDispPrior <- median(od_expressed[od_expressed > 1])
    if (!is.finite(OverDispPrior) || OverDispPrior < 1) OverDispPrior <- 1
    DFPrior <- 10
    OverDisp[expressed] <- (DFPrior * OverDispPrior + DF[expressed] * od_expressed) /
      (DFPrior + DF[expressed])
    OverDisp[expressed] <- pmax(OverDisp[expressed], 1)
  }
  
  .log("[OD Bootstrap] range: min=%.3f, med=%.3f, max=%.3f, pct>1.5=%.1f%%",
       min(OverDisp), median(OverDisp), max(OverDisp), mean(OverDisp > 1.5) * 100)
  
  list(OverDisp = as.numeric(OverDisp), DF = as.integer(DF),
       features = genes, expr_cells = as.numeric(expr_cells))
}


#' Estimate overdispersion using method-of-moments
#'
#' Computes per-gene overdispersion estimates from a counts matrix
#' using a method-of-moments approach. Optionally normalizes counts
#' by size factors before estimation.
#'
#' @param counts A numeric or sparse matrix of gene-by-cell counts.
#' @param min_cells Minimum number of cells with nonzero expression
#'   required for a gene to be included in estimation. Transcripts with fewer
#'   than this number are assigned an overdispersion of 1.
#' @param size_factors Optional numeric vector of size factors for each
#'   cell. If \code{NULL}, they are computed as the column sums of
#'   \code{counts}, median-scaled to 1.
#'
#' @details
#' For each gene, nonzero expression values are used to compute the mean
#' and variance. A pseudocount of 0.1 is added to both before computing
#' overdispersion (\eqn{OD = (var + 0.1) / (mean + 0.1)}). The estimate
#' is then shrunk toward 1 based on the number of expressing cells, with
#' stronger shrinkage for small sample sizes.
#'
#' @return A numeric vector of overdispersion estimates (length equal to
#'   the number of genes).
#'
#' @importFrom Matrix colSums Diagonal
#' @export


estimate_overdispersion_moments <- function(counts, min_cells = 10, size_factors = NULL) {
  if (is.null(size_factors)) {
    size_factors <- Matrix::colSums(counts)
    size_factors <- size_factors / median(size_factors[size_factors > 0])
  }
  counts_norm <- counts %*% Matrix::Diagonal(x = 1 / size_factors)
  
  n_genes <- nrow(counts_norm)
  od_estimates <- numeric(n_genes)
  
  for (g in 1:n_genes) {
    gene_expr <- as.numeric(counts_norm[g, ])
    expr_nonzero <- gene_expr[gene_expr > 0]
    if (length(expr_nonzero) < min_cells) {
      od_estimates[g] <- 1
    } else {
      mu <- mean(expr_nonzero)
      var_obs <- var(expr_nonzero)
      od <- (var_obs + 0.1) / (mu + 0.1)
      n <- length(expr_nonzero)
      confidence <- n / (n + 30)
      od_shrunk <- 1 + confidence * (od - 1)
      od_estimates[g] <- max(od_shrunk, 1)
    }
  }
  
  .log("[OD Moments] range: min=%.3f, med=%.3f, max=%.3f, pct>1.5=%.1f%%",
       min(od_estimates), median(od_estimates), max(od_estimates),
       mean(od_estimates > 1.5) * 100)
  od_estimates
}


#' Calculate complexity-based prior for overdispersion
#'
#' Computes a per-gene prior scaling factor for overdispersion estimates
#' based on transcript complexity, using either a transcript info table
#' or a GTF annotation file.
#'
#' @param gene_names Character vector of gene IDs for which to compute the prior.
#' @param gtf_file Optional path to a GTF file containing transcript annotations.
#'   If provided and \code{transcript_info} is not given, transcript counts per
#'   gene will be computed from this file.
#' @param transcript_info Optional path to a tab-delimited file containing
#'   transcript metadata with a \code{gene_id} column. Used to count transcripts
#'   per gene.
#'
#' @details
#' For each gene, the number of transcripts is counted, and the prior scaling
#' factor is computed as:
#' \deqn{prior = 1 + 0.3 \times \log1p(n\_transcripts - 1)}
#' Genes absent from the annotation are assigned a prior of 1.0.
#'
#' @return A numeric vector of prior scaling factors (length equal to the number
#'   of genes), where larger values correspond to higher transcript complexity.
#'
#' @importFrom rtracklayer import
#' @export


calculate_complexity_prior <- function(gene_names, gtf_file = NULL, transcript_info = NULL) {
  n_genes <- length(gene_names)
  prior <- rep(1.0, n_genes)
  
  if (!is.null(transcript_info) && file.exists(transcript_info)) {
    tx_info <- read.table(transcript_info, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    tx_per_gene <- table(tx_info$gene_id)
    for (i in 1:n_genes) {
      gene <- gene_names[i]
      if (gene %in% names(tx_per_gene)) {
        n_transcripts <- as.numeric(tx_per_gene[gene])
        prior[i] <- 1 + 0.3 * log1p(n_transcripts - 1)
      }
    }
  } else if (!is.null(gtf_file) && file.exists(gtf_file)) {
    gtf <- rtracklayer::import(gtf_file)
    transcripts <- gtf[gtf$type == "transcript"]
    tx_per_gene <- table(transcripts$gene_id)
    for (i in 1:n_genes) {
      gene <- gene_names[i]
      if (gene %in% names(tx_per_gene)) {
        n_transcripts <- as.numeric(tx_per_gene[gene])
        prior[i] <- 1 + 0.3 * log1p(n_transcripts - 1)
      }
    }
  }
  
  .log("[OD Prior] range: min=%.3f, med=%.3f, max=%.3f, pct>1.5=%.1f%%",
       min(prior), median(prior), max(prior), mean(prior > 1.5) * 100)
  prior
}


#' Integrated overdispersion estimation
#'
#' Estimates per-gene overdispersion by integrating three complementary methods:
#' (1) bootstrap-based estimation from Alevin output, (2) method-of-moments
#' estimation from normalized counts, and (3) a complexity-based prior
#' from transcript annotation. Produces a weighted geometric mean of the three
#' methods, followed by empirical Bayes shrinkage.
#'
#' @param counts A sparse or dense numeric matrix of raw counts (genes x cells).
#' @param alevin_dir Path to the directory containing Alevin output with
#'   bootstrap replicates.
#' @param n_boot Integer, number of bootstrap replicates per cell in the
#'   Alevin output.
#' @param gene_names Optional character vector of gene names (length equal to
#'   number of rows in \code{counts}). If \code{NULL}, rownames of \code{counts}
#'   are used.
#' @param gtf_file Optional path to a GTF file for computing transcript complexity.
#' @param transcript_info Optional path to a transcript info table for computing
#'   transcript complexity (must contain \code{gene_id} column).
#' @param method_weights Named numeric vector of weights for combining the three
#'   methods (\code{"bootstrap"}, \code{"moments"}, \code{"prior"}). Will be
#'   normalized to sum to 1.
#' @param min_cells_expr Minimum number of expressing cells required for a gene
#'   to be included in moment-based and bootstrap-based estimates.
#' @param n_cores Number of parallel cores to use for bootstrap-based estimation.
#' @param debug_first_block Logical, if \code{TRUE} only the first bootstrap block
#'   is processed (for debugging).
#'
#' @details
#' The function integrates multiple sources of information on technical noise:
#' \enumerate{
#'   \item \strong{Bootstrap-based}: computes overdispersion from per-cell bootstrap
#'     replicates in Alevin output.
#'   \item \strong{Moment-based}: uses nonzero normalized counts to estimate variance
#'     inflation relative to the mean.
#'   \item \strong{Complexity prior}: uses transcript-level annotation to adjust
#'     expected variance inflation based on gene complexity.
#' }
#' The weighted geometric mean of these three estimates is calculated, then shrunk
#' toward the robust median for well-observed genes (\eqn{\ge 50} cells).
#'
#' @return A list with elements:
#' \describe{
#'   \item{\code{overdispersion}}{Final integrated overdispersion estimates.}
#'   \item{\code{od_bootstrap}}{Bootstrap-based estimates.}
#'   \item{\code{od_moments}}{Moment-based estimates.}
#'   \item{\code{od_prior}}{Complexity prior values.}
#'   \item{\code{expr_cells}}{Number of cells expressing each gene.}
#'   \item{\code{features}}{Gene names.}
#' }
#'
#' @importFrom Matrix colSums Diagonal
#' @export


estimate_overdispersion_integrated <- function(
    counts,
    alevin_dir,
    n_boot,
    gene_names     = NULL,
    gtf_file       = NULL,
    transcript_info = NULL,
    method_weights = c(bootstrap = 0.3, moments = 0.5, prior = 0.2),
    min_cells_expr = 10,
    n_cores        = 4,
    debug_first_block = FALSE
) {
  if (is.null(gene_names)) gene_names <- rownames(counts)
  n_genes <- length(gene_names)
  .log("Starting integrated overdispersion estimation for %d genes", n_genes)
  
  # Method 1: Bootstrap-based (sparse-aware)
  od_boot_result <- compute_overdisp_sparse_aware(
    alevin_dir = alevin_dir,
    n_boot = n_boot,
    min_cells_expr = min_cells_expr,
    pseudocount = 0.1,
    n_cores = n_cores,
    verbose = TRUE,
    debug_first_block = debug_first_block
  )
  od_bootstrap <- od_boot_result$OverDisp
  
  # Method 2: Moment-based
  od_moments <- estimate_overdispersion_moments(
    counts = counts,
    min_cells = min_cells_expr
  )
  
  # Method 3: Complexity prior
  od_prior <- calculate_complexity_prior(
    gene_names = gene_names,
    gtf_file = gtf_file,
    transcript_info = transcript_info
  )
  
  # Normalize weights
  method_weights <- method_weights / sum(method_weights)
  
  # Weighted geometric mean (ratio-scale)
  log_od_integrated <- method_weights["bootstrap"] * log(pmax(od_bootstrap, 1)) +
    method_weights["moments"]   * log(pmax(od_moments, 1))   +
    method_weights["prior"]     * log(pmax(od_prior, 1))
  od_integrated <- exp(log_od_integrated)
  
  # Final EB shrinkage toward robust center for sufficiently observed genes
  expr_cells <- od_boot_result$expr_cells
  high_conf <- which(expr_cells >= 50)
  if (length(high_conf) > 0) {
    od_median_robust <- median(od_integrated[high_conf])
    for (i in seq_len(n_genes)) {
      shrink_weight <- min(expr_cells[i] / 100, 1)  # less shrinkage with more data
      od_integrated[i] <- shrink_weight * od_integrated[i] +
        (1 - shrink_weight) * od_median_robust
    }
  }
  
  # Stability clamps
  od_integrated <- pmax(od_integrated, 1)
  od_integrated <- pmin(od_integrated, 100)
  
  .log("[OD Integrated] Final range: min=%.3f, med=%.3f, max=%.3f",
       min(od_integrated), median(od_integrated), max(od_integrated))
  .log("[OD Integrated] Transcripts with OD > 1.5: %.1f%%, > 2: %.1f%%",
       mean(od_integrated > 1.5) * 100, mean(od_integrated > 2) * 100)
  
  list(
    overdispersion = od_integrated,
    od_bootstrap   = od_bootstrap,
    od_moments     = od_moments,
    od_prior       = od_prior,
    expr_cells     = expr_cells,
    features       = gene_names
  )
}


#' Read sample count data and estimate integrated overdispersion
#'
#' Loads a gene-by-cell count matrix for a given sample from Alevin output,
#' then computes integrated overdispersion estimates by combining bootstrap,
#' moment-based, and complexity prior methods.
#'
#' @param sample_id Character, sample identifier (subdirectory name under
#'   \code{base_dir}).
#' @param base_dir Character, path to the base directory containing sample
#'   subdirectories with Alevin output.
#' @param n_boot Integer, number of bootstrap replicates per cell in the
#'   Alevin output.
#' @param gtf_file Optional path to a GTF file for transcript complexity
#'   calculation.
#' @param transcript_info Optional path to a transcript info table for
#'   complexity calculation (must contain \code{gene_id} column).
#' @param block_cells Integer, number of cells per processing block for
#'   bootstrap-based estimation.
#' @param n_cores Integer, number of parallel cores to use for processing.
#'   Defaults to one less than the number of available cores.
#' @param omp_threads Integer, number of OpenMP threads to use for C++ code
#'   (if applicable).
#' @param debug_first_block Logical, if \code{TRUE} only the first bootstrap
#'   block is processed (for debugging).
#'
#' @details
#' This function reads the sparse counts matrix (\code{quants_mat.gz}) along
#' with its associated feature and cell name files from the specified Alevin
#' directory. The counts are loaded via \code{\link{read_eds_gc}} to ensure
#' dimension checks and logging.
#'
#' After loading, \code{\link{estimate_overdispersion_integrated}} is called
#' to compute integrated overdispersion estimates, which combine:
#' \enumerate{
#'   \item Bootstrap-based estimation from per-cell Alevin replicates.
#'   \item Method-of-moments estimation from normalized counts.
#'   \item A complexity prior from transcript annotation.
#' }
#'
#' @return A list with elements:
#' \describe{
#'   \item{\code{counts}}{Sparse counts matrix (genes x cells).}
#'   \item{\code{od}}{Result list from \code{\link{estimate_overdispersion_integrated}}.}
#'   \item{\code{feats}}{Character vector of gene names.}
#'   \item{\code{cells}}{Character vector of cell barcodes.}
#' }
#'
#' @export


read_sample_data_improved <- function(
    sample_id,
    base_dir,
    n_boot,
    gtf_file = NULL,
    transcript_info = NULL,
    block_cells = 128L,
    n_cores     = max(1L, parallel::detectCores() - 1L),
    omp_threads = 1L,
    debug_first_block = FALSE
) {
  t0 <- Sys.time()
  alevin_dir <- file.path(base_dir, sample_id, "alevin")
  .log("[%s] Starting in %s", sample_id, alevin_dir)
  
  # Read count matrix
  f_counts <- file.path(alevin_dir, "quants_mat.gz")
  f_cols   <- file.path(alevin_dir, "quants_mat_cols.txt")
  f_rows   <- file.path(alevin_dir, "quants_mat_rows.txt")
  stopifnot(file.exists(f_counts), file.exists(f_cols), file.exists(f_rows))
  
  feats <- scan(f_cols, what = "character", quiet = TRUE)
  cells <- scan(f_rows, what = "character", quiet = TRUE)
  .log("[%s] Features=%d (genes) Cells=%d", sample_id, length(feats), length(cells))
  
  counts <- read_eds_gc(
    f_counts,
    n_feat  = length(feats),
    n_cells = length(cells),
    label   = sprintf("%s:counts", sample_id)
  )
  stopifnot(nrow(counts) == length(feats), ncol(counts) == length(cells))
  rownames(counts) <- feats
  colnames(counts) <- paste0(cells, "_", sample_id)
  
  # Integrated overdispersion estimation
  od_result <- estimate_overdispersion_integrated(
    counts = counts,
    alevin_dir = alevin_dir,
    n_boot = n_boot,
    gene_names = feats,
    gtf_file = gtf_file,
    transcript_info = transcript_info,
    method_weights = c(bootstrap = 0.3, moments = 0.5, prior = 0.2),
    min_cells_expr = 10,
    n_cores = n_cores,
    debug_first_block = debug_first_block
  )
  
  .log("[%s] Data loading and integrated OD calculation complete", sample_id)
  
  list(
    counts = counts,
    od     = od_result,
    feats  = feats,
    cells  = cells
  )
}


#' Create a Seurat Object with a Corrected Assay
#'
#' Generates a Seurat object containing both the original raw counts and an
#' overdispersion-corrected count assay. Default correction is based on bootstrap
#' estimates from \code{\link{compute_overdisp_sparse_aware}}.
#'
#' @param sample_id Character string; a unique identifier for the sample, used
#'   for labeling the Seurat object project.
#' @param counts A sparse matrix of raw gene-by-cell counts.
#' @param od A list returned by \code{\link{estimate_overdispersion_integrated}},
#'   containing the final integrated `overdispersion` estimates and its
#'   component values.
#' @param feats Character vector of feature (gene) names, which must match the
#'   rownames of the `counts` matrix.
#' @param cells Character vector of cell barcodes, which must match the
#'   colnames of the `counts` matrix.
#'
#' @details
#' This function creates a Seurat object with the raw counts stored in the default
#' `RNA` assay. It then calculates corrected counts by dividing the raw counts
#' for each gene by the corresponding bootstrap overdispersion factor from the
#' `od` object. These corrected counts are stored in a new assay named `RNA_corr`.
#'
#' Both assays are populated with rich feature-level metadata, including the
#' integrated overdispersion value and its components (bootstrap, moments, prior),
#' which can be accessed via `Misc(seu[["assay"]], "feature_meta")`.
#'
#' @return A Seurat object with two assays:
#' \describe{
#'   \item{`RNA`}{The original, uncorrected count data.}
#'   \item{`RNA_corr`}{The overdispersion-corrected count data. This is set as the default assay.}
#' }
#' The object's `misc` slot also stores metadata about the correction process.
#'
#' @importFrom Seurat CreateSeuratObject CreateAssayObject DefaultAssay
#' @importFrom SeuratObject Misc
#' @export


process_and_create_seurat_corrected_improved <- function(
    sample_id,
    counts,
    od,
    feats,
    cells
) {
  t0 <- Sys.time()
  .log("[%s] Starting Seurat object creation with IMPROVED scaling", sample_id)
  
  scaling_factors <- od$od_bootstrap
  if (length(scaling_factors) != nrow(counts)) {
    stop("Mismatch between overdispersion length and count matrix rows")
  }
  
  .log("[%s] Integrated OD range: min=%.3f, median=%.3f, max=%.3f",
       sample_id, min(scaling_factors), median(scaling_factors), max(scaling_factors))
  .log("[%s] Transcripts with OD > 1.5: %.1f%%, > 2: %.1f%%",
       sample_id, mean(scaling_factors > 1.5) * 100, mean(scaling_factors > 2) * 100)
  
  inv_scaling <- 1 / scaling_factors
  counts_corr <- Matrix::Diagonal(x = inv_scaling) %*% counts
  rownames(counts_corr) <- rownames(counts)
  colnames(counts_corr) <- colnames(counts)
  
  seu <- Seurat::CreateSeuratObject(counts = counts, project = sample_id, assay = "RNA")
  seu[["RNA_corr"]] <- Seurat::CreateAssayObject(counts = counts_corr)
  
  fm <- data.frame(
    OverDisp_integrated = od$overdispersion,
    OverDisp_bootstrap  = od$od_bootstrap,
    OverDisp_moments    = od$od_moments,
    OverDisp_prior      = od$od_prior,
    expressing_cells    = od$expr_cells,
    scaling_factor      = scaling_factors,
    inv_scaling         = inv_scaling,
    row.names = feats
  )
  SeuratObject::Misc(seu[["RNA"]], "feature_meta")      <- fm
  SeuratObject::Misc(seu[["RNA_corr"]], "feature_meta") <- fm
  SeuratObject::Misc(seu, "od_method") <- "bootstrap"
  
  DefaultAssay(seu) <- "RNA_corr"
  seu$sample <- sample_id
  
  .log("[%s] Done in %.1fs | %d genes x %d cells", sample_id,
       as.numeric(difftime(Sys.time(), t0, "secs")), nrow(seu), ncol(seu))
  seu
}


#' Harmonize Feature Metadata in Seurat Objects
#'
#' Aligns feature metadata across multiple assays in a list of Seurat objects.
#' This function is useful for replacing existing feature IDs (e.g., Ensembl IDs)
#' with more human-readable ones (e.g., gene symbols) from an external file, while
#' preserving the original IDs in a separate column.
#'
#' @param objs A list of Seurat objects to be modified.
#' @param file_path Path to a tab-separated file where the first column contains
#'   the new feature IDs (e.g., gene symbols). The order of IDs in this file must
#'   match the order of features in the Seurat objects.
#' @param assays A character vector specifying which assays within each Seurat
#'   object to update (e.g., `c("RNA", "RNA_corr")`).
#'
#' @details
#' For each specified assay in each Seurat object, this function accesses the
#' feature metadata stored in `object[[assay]]@misc$feature_meta`. It performs
#' the following actions:
#' \enumerate{
#'   \item Replaces the `rownames` of the `feature_meta` data frame with the IDs
#'         from the supplied file.
#'   \item Adds a new column named `feature_id` to the `feature_meta` data frame,
#'         containing the original feature names from the assay's count matrix.
#'   \item Clears the `rownames` of the modified `feature_meta` data frame.
#' }
#' This ensures that the original feature identifiers are not lost and can be
#' referenced later.
#'
#' @return The input list of Seurat objects, with the `feature_meta` in the
#'   specified assays modified in place.
#' @export


set_feature_metadata <- function(objs, file_path, assays = c("RNA", "RNA_corr")) {
  # Read feature IDs from the provided file path once
  feature_ids <- as.character(read.table(file_path, sep="\t", header=T)[,1])
  # feature_ids <- unique(as.character(read.table(file_path, sep="\t", header=F)[,2]))
  stopifnot(is.list(objs))
  
  # Process each Seurat object in the list
  lapply(objs, function(o) {
    # Process each specified assay within the object
    for (assay in assays) {
      # Check if assay exists
      if (!assay %in% names(o@assays)) {
        warning(sprintf("[set_feature_metadata] Assay '%s' not found; skipping.", assay))
        next
      }
      
      # Retrieve feature_meta and check for null
      fm <- SeuratObject::Misc(o[[assay]], "feature_meta")
      if (is.null(fm)) {
        warning(sprintf("[set_feature_metadata] feature_meta is NULL in assay '%s'; skipping.", assay))
        next
      }
      fm <- if (is.data.frame(fm)) fm else as.data.frame(fm, stringsAsFactors = FALSE)
      
      # Check for length mismatch with file-based feature IDs
      if (length(feature_ids) != nrow(fm)) {
        stop(sprintf("[set_feature_metadata] Length mismatch in assay '%s': length(feature_ids)=%d, nrow(feature_meta)=%d",
                     assay, length(feature_ids), nrow(fm)))
      }
      
      # Assign file-based IDs to rownames of feature_meta
      rownames(fm) <- make.unique(feature_ids)
      
      # Get the feature names from the assay's counts matrix
      assay_features <- rownames(GetAssayData(o, assay = assay, layer = "counts"))
      
      # Check for length mismatch with assay features
      if (nrow(fm) != length(assay_features)) {
        warning(sprintf("[set_feature_metadata] '%s': nrow(feature_meta)=%d != nfeatures=%d; skipping.",
                        assay, nrow(fm), length(assay_features)))
        next
      }
      
      # Add a new column with the assay's feature names
      fm$feature_id <- assay_features
      
      # Clear the rownames as they are no longer the primary identifier
      rownames(fm) <- NULL
      
      # Write the modified feature_meta back to the Seurat object
      SeuratObject::Misc(o[[assay]], "feature_meta") <- fm
    }
    # Return the modified Seurat object
    o
  })
}


#' Timestamped logging utility
#'
#' Prints a formatted log message to the console, prefixed with the current
#' time in \code{HH:MM:SS} format.
#'
#' @param fmt Character string. Message template, optionally containing
#'   \code{\%} format specifiers for use with \code{\link[base]{sprintf}}.
#' @param ... Optional arguments to be inserted into \code{fmt} via
#'   \code{\link[base]{sprintf}}, or appended directly if no format specifiers
#'   are present.
#'
#' @details
#' If \code{fmt} contains \code{\%} format specifiers, the function uses
#' \code{\link[base]{sprintf}} to substitute in the provided arguments.
#' Otherwise, all arguments are concatenated into a single message string.
#' The output is automatically prepended with a timestamp and printed
#' to the console with \code{\link[base]{message}}.
#'
#' @return Invisibly returns \code{NULL}. Called for its side effect of
#'   printing a log message.
#'
#' @export

.log <- function(fmt, ...) {
  args <- list(...)
  msg <-
    if (length(args) && grepl("%", fmt, fixed = TRUE)) do.call(sprintf, c(fmt, args))
  else paste0(c(fmt, unlist(args)), collapse = "")
  message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), msg))
}


#' Extract a Named Feature Vector from Seurat Objects
#'
#' Searches through a list of Seurat objects to find and extract a specific
#' metadata vector from the feature metadata stored in one of the assays.
#'
#' @param seurat_list A list of Seurat objects.
#' @param vector_name The character name of the column to extract from the
#'   feature metadata data frame (e.g., "scaling_factor").
#' @param prefer_assays A character vector of assay names to search, in order
#'   of preference (e.g., `c("RNA_corr", "RNA")`).
#'
#' @details
#' The function iterates through each object and then through the preferred assays.
#' It stops and returns the vector as soon as it is found in the first matching
#' location (`object@assays[[assay]]@misc$feature_meta`).
#'
#' @return A named numeric or character vector, where names are the feature IDs.
#'   Returns `NULL` and a warning if the vector is not found in any of the
#'   specified locations.
#'
#' @seealso \code{\link{calculate_within_sample_bcv}}, \code{\link{calculate_bcv}}
#' @export


extract_feature_vector <- function(seurat_list, vector_name,
                                   prefer_assays = c("RNA_corr","RNA")) {
  .log("Searching for feature vector '%s' across %d objects", vector_name, length(seurat_list))
  
  for (i in seq_along(seurat_list)) {
    obj <- seurat_list[[i]]
    for (assay_nm in prefer_assays) {
      if (assay_nm %in% names(obj@assays)) {
        fm <- obj@assays[[assay_nm]]@misc[["feature_meta"]]
        if (!is.null(fm) && vector_name %in% names(fm)) {
          vec <- fm[[vector_name]]
          names(vec) <- fm[["feature_id"]]
          .log("Found '%s' in object %d assay %s for %d genes", vector_name, i, assay_nm, length(vec))
          return(vec)
        }
      }
    }
  }
  warning(sprintf("Vector '%s' not found in feature_meta of any object/assay", vector_name))
  NULL
}


#' Extract Overdispersion Vector
#'
#' A wrapper function to extract the default overdispersion vector
#' ("scaling_factor") from a list of Seurat objects.
#'
#' @param seurat_list A list of Seurat objects.
#'
#' @return A named numeric vector of overdispersion values. Returns `NULL` if
#'   "scaling_factor" is not found in the feature metadata.
#'
#' @details This function is a convenient helper that uses
#'   \code{\link{extract_feature_vector}} to retrieve overdispersion estimates.
#'
#' @export
extract_overdispersion <- function(seurat_list) {
  v <- extract_feature_vector(seurat_list, "scaling_factor")
  if (!is.null(v)) return(v)
}


#' Subsample cells into groups for BCV calculation
#'
#' Randomly partitions cells from a Seurat object into a specified number of
#' groups for biological coefficient of variation (BCV) estimation.
#'
#' @param seurat_obj A \code{\link[Seurat]{Seurat}} object containing single-cell
#'   expression data.
#' @param n_groups Integer. Desired number of cell groups to create
#'   (default: 10).
#' @param cells_per_group Integer. Number of cells per group. If \code{NULL}
#'   (default), this is calculated as \code{max(min_cells_per_group, floor(n_cells / n_groups))}.
#' @param min_cells_per_group Minimum number of cells per group if
#'   \code{cells_per_group} is not specified (default: 50).
#' @param seed Integer random seed for reproducibility (default: 42).
#'
#' @details
#' The function shuffles cell IDs and splits them into approximately equal-sized
#' groups. If the requested \code{n_groups} would result in fewer cells per group
#' than \code{min_cells_per_group}, it is adjusted automatically.
#'
#' This grouping is intended for downstream BCV (biological coefficient of
#' variation) calculations, ensuring each group has enough cells for reliable
#' dispersion estimation.
#'
#' @return A named list where each element contains the cell IDs for a group.
#'
#' @export


subsample_cells_for_bcv <- function(seurat_obj, 
                                    n_groups = 10,
                                    cells_per_group = NULL,
                                    min_cells_per_group = 50,
                                    seed = 42) {
  
  set.seed(seed)
  n_cells <- ncol(seurat_obj)
  
  if (is.null(cells_per_group)) {
    cells_per_group <- max(min_cells_per_group, floor(n_cells / n_groups))
  }
  
  max_groups <- floor(n_cells / cells_per_group)
  if (n_groups > max_groups) {
    n_groups <- max_groups
    .log("Adjusted n_groups to %d based on available cells", n_groups)
  }
  if (n_groups < 2) stop("Not enough cells to create at least 2 groups for BCV calculation")
  
  cell_ids <- colnames(seurat_obj)
  shuffled <- sample(cell_ids)
  
  groups <- vector("list", n_groups)
  for (i in seq_len(n_groups)) {
    start_idx <- (i - 1) * cells_per_group + 1
    end_idx <- min(i * cells_per_group, length(shuffled))
    groups[[paste0("Group_", i)]] <- shuffled[start_idx:end_idx]
  }
  
  .log("Created %d groups with ~%d cells each", n_groups, cells_per_group)
  groups
}


#' Create within-sample pseudo-bulk counts
#'
#' Aggregates single-cell counts from specified cell groups into a
#' pseudo-bulk count matrix for a given assay.
#'
#' @param seurat_obj A \code{\link[Seurat]{Seurat}} object containing single-cell
#'   expression data.
#' @param cell_groups A named list of character vectors, where each vector
#'   contains the cell barcodes for one group (e.g., from
#'   \code{\link{subsample_cells_for_bcv}}).
#' @param assay_name Character scalar. The assay to use for counts extraction
#'   (default: \code{"RNA"}).
#'
#' @details
#' This function takes a set of predefined groups of cells from a Seurat object
#' and sums their counts to produce a pseudo-bulk matrix. This can be used for
#' downstream analysis such as biological coefficient of variation (BCV)
#' estimation or differential expression at the group level.
#'
#' @return A sparse matrix (\code{dgCMatrix}) with genes as rows and groups as
#'   columns, containing summed counts for each group.
#'
#' @importFrom Seurat GetAssayData
#' @importFrom Matrix rowSums
#'
#' @export


create_within_sample_pseudobulk <- function(seurat_obj, cell_groups, assay_name = "RNA") {
  counts <- GetAssayData(seurat_obj, assay = assay_name, layer = "counts")
  pb_matrix <- matrix(0, nrow = nrow(counts), ncol = length(cell_groups))
  rownames(pb_matrix) <- rownames(counts)
  colnames(pb_matrix) <- names(cell_groups)
  for (i in seq_along(cell_groups)) {
    group_cells <- intersect(cell_groups[[i]], colnames(counts))
    if (length(group_cells) > 0) {
      pb_matrix[, i] <- Matrix::rowSums(counts[, group_cells, drop = FALSE])
    }
  }
  as(pb_matrix, "dgCMatrix")
}


#' Calculate within-sample BCV using pseudo-bulk groups
#'
#' Forms within-sample pseudo-bulk groups from a Seurat object and estimates
#' the biological coefficient of variation (BCV) for each gene using edgeR.
#'
#' @param seurat_obj A \code{\link[Seurat]{Seurat}} object containing single-cell
#'   expression data.
#' @param assay_name Character scalar. The assay to use for counts extraction
#'   (default: \code{"RNA"}).
#' @param n_groups Integer. Number of cell groups to form within the sample
#'   (default: 10).
#' @param cells_per_group Integer. Number of cells per group; if \code{NULL},
#'   will be determined based on \code{n_groups} and available cells.
#' @param min_cells_per_group Integer. Minimum number of cells per group
#'   (default: 50).
#' @param min_counts Integer. Minimum total counts per gene across all groups
#'   for the gene to be retained (default: 10).
#' @param robust Logical. Whether to use robust dispersion estimation in
#'   \code{\link[edgeR]{estimateDisp}} (default: TRUE).
#' @param seed Integer. Random seed for reproducibility (default: 42).
#' @param feature_vector Optional numeric vector of per-gene values (named by
#'   gene) to be applied as an offset in the model (e.g., overdispersion or
#'   scaling ratios).
#' @param vector_as Character. Interpretation of \code{feature_vector}: either
#'   \code{"od"} for overdispersion or \code{"ratio"} for scaling factors.
#'
#' @details
#' This function:
#' \enumerate{
#'   \item Subsamples cells into \code{n_groups} groups (via
#'     \code{\link{subsample_cells_for_bcv}}).
#'   \item Aggregates counts into group-level pseudo-bulk profiles (via
#'     \code{\link{create_within_sample_pseudobulk}}).
#'   \item Filters genes by total counts across groups.
#'   \item Optionally applies a per-gene offset from a supplied feature vector.
#'   \item Estimates tagwise, common, and trended dispersions using edgeR.
#' }
#'
#' The biological coefficient of variation (BCV) is the square root of the
#' dispersion estimate and can be used to assess variability independent of
#' mean expression.
#'
#' @return A list with:
#' \describe{
#'   \item{genes}{Character vector of retained gene IDs.}
#'   \item{groups}{Character vector of group names.}
#'   \item{n_cells_total}{Total number of cells in the Seurat object.}
#'   \item{n_groups}{Number of groups after filtering.}
#'   \item{cells_per_group}{Number of cells per group.}
#'   \item{ave_log_cpm}{Average log2 counts per million (CPM) per gene.}
#'   \item{bcv_tagwise}{Tagwise BCV estimates.}
#'   \item{bcv_common}{Common BCV estimate.}
#'   \item{bcv_trend}{Trended BCV estimates.}
#'   \item{n_genes}{Number of genes retained.}
#' }
#'
#' @importFrom Matrix rowSums colSums
#' @importFrom edgeR DGEList calcNormFactors estimateDisp
#' @importFrom stats model.matrix
#'
#' @export


calculate_within_sample_bcv <- function(seurat_obj,
                                        assay_name = "RNA",
                                        n_groups = 10,
                                        cells_per_group = NULL,
                                        min_cells_per_group = 50,
                                        min_counts = 10,
                                        robust = TRUE,
                                        seed = 42,
                                        feature_vector = NULL,
                                        vector_as = c("od","ratio")) {
  vector_as <- match.arg(vector_as)
  
  # Form within-sample groups and pseudobulk
  cell_groups <- subsample_cells_for_bcv(
    seurat_obj, n_groups, cells_per_group, min_cells_per_group, seed
  )
  pb <- create_within_sample_pseudobulk(seurat_obj, cell_groups, assay_name)
  
  # Gene filtering
  keep <- Matrix::rowSums(pb) >= min_counts
  pb_filtered <- pb[keep, , drop = FALSE]
  if (nrow(pb_filtered) < 100) {
    warning("Very few genes passed filtering. Consider reducing min_counts.")
  }
  
  # Ensure unique group names to avoid "Repeated column names" warnings
  colnames(pb_filtered) <- make.unique(colnames(pb_filtered))
  
  # Drop groups with zero library size after filtering (fix: use elementwise '&', not '&&')
  lib_sizes <- Matrix::colSums(pb_filtered)
  keep_cols <- is.finite(lib_sizes) & (lib_sizes > 0)
  if (any(!keep_cols)) {
    .log("Dropping %d/%d empty groups (zero library after filtering)", sum(!keep_cols), length(keep_cols))
    pb_filtered <- pb_filtered[, keep_cols, drop = FALSE]
  }
  if (ncol(pb_filtered) < 2) {
    stop("After filtering, fewer than 2 non-empty groups remain; try lowering min_counts or cells_per_group.")
  }
  
  # edgeR objects
  y <- edgeR::DGEList(counts = pb_filtered)
  y <- edgeR::calcNormFactors(y, method = "TMM")
  design <- stats::model.matrix(~ 1, data = data.frame(group = colnames(pb_filtered)))
  
  # Optional offsets from a per-gene feature vector
  if (!is.null(feature_vector)) {
    fv <- feature_vector[rownames(y$counts)]
    if (any(is.na(fv))) {
      .log("Some genes missing '%s' vector values; setting to neutral (1.0)", vector_as)
      fv[is.na(fv)] <- 1
    }
    # Base offset
    base_off <- matrix(
      log(y$samples$lib.size * y$samples$norm.factors),
      nrow = nrow(y$counts), ncol = ncol(y$counts), byrow = TRUE
    )
    # Add per-gene offset
    fv <- pmax(fv, 1e-12)
    add_off <- matrix(log(fv), nrow = length(fv), ncol = ncol(y$counts), byrow = FALSE)
    y$offset <- base_off + add_off
    storage.mode(y$offset) <- "double"
    .log("Applied within-sample offsets using feature vector as '%s'", vector_as)
  }
  
  # Dispersion/BCV
  y <- edgeR::estimateDisp(y, design = design, robust = robust)
  
  list(
    genes = rownames(y$counts),
    groups = colnames(y$counts),
    n_cells_total = ncol(seurat_obj),
    n_groups = ncol(y$counts),
    cells_per_group = sapply(cell_groups, length),
    ave_log_cpm = y$AveLogCPM,
    bcv_tagwise = sqrt(y$tagwise.dispersion),
    bcv_common = sqrt(y$common.dispersion),
    bcv_trend = sqrt(y$trended.dispersion),
    n_genes = nrow(y$counts)
  )
}


#' Analyze Within-Sample BCV Across Multiple Seurat Objects
#'
#' This function computes within-sample biological coefficient of variation (BCV) 
#' metrics for each Seurat object in a list, optionally across multiple assays 
#' (e.g., raw and overdispersion-corrected counts). BCV is estimated via 
#' \code{\link{calculate_within_sample_bcv}} using pseudobulk replicates formed 
#' by random subsampling of cells within each sample.
#'
#' @param seurat_list A named or unnamed list of Seurat objects to process.
#' @param assay_names Character vector of assay names to evaluate (default: \code{c("RNA","RNA_corr")}).
#' @param n_groups Integer, number of pseudobulk groups to create per sample (default: 10).
#' @param cells_per_group Optional integer specifying the number of cells per group; 
#'   if \code{NULL}, determined from \code{n_groups} and total cells.
#' @param min_cells_per_group Minimum number of cells allowed per pseudobulk group (default: 50).
#' @param min_counts Minimum total counts per gene across all groups to retain it for BCV 
#'   calculation (default: 10).
#' @param robust Logical, whether to use robust dispersion estimation in \code{edgeR} (default: TRUE).
#' @param seed Integer seed for reproducible grouping (default: 42).
#' @param offset_vector_name Optional character string naming a feature metadata vector 
#'   to use as an offset in BCV estimation (e.g., "scaling_factor").
#' @param offset_vector_as Character scalar, either \code{"od"} (interpret offset vector 
#'   as overdispersion factors) or \code{"ratio"} (interpret as multiplicative ratios 
#'   to adjust counts). Default is \code{"od"}.
#'
#' @return A named list, with one element per assay. Each element is itself a list of BCV 
#'   results for each sample, as returned by \code{\link{calculate_within_sample_bcv}}, 
#'   augmented with the \code{sample_name}.
#'
#' @details
#' For each assay in each sample, cells are randomly assigned to \code{n_groups} pseudobulk 
#' replicates. BCV is calculated across these replicates using edgeR's dispersion estimation.
#' If \code{offset_vector_name} is provided, the function will attempt to extract the 
#' corresponding per-gene vector from the feature metadata of the first matching assay 
#' found in the Seurat objects, and will pass it to \code{\link{calculate_within_sample_bcv}}.
#'
#' @examples
#' \dontrun{
#' results <- analyze_within_sample_bcv(
#'   seurat_list = list(sample1 = seu1, sample2 = seu2),
#'   assay_names = c("RNA", "RNA_corr"),
#'   n_groups = 8,
#'   offset_vector_name = "scaling_factor"
#' )
#' }
#'
#' @export


analyze_within_sample_bcv <- function(seurat_list,
                                      assay_names = c("RNA", "RNA_corr"),
                                      n_groups = 10,
                                      cells_per_group = NULL,
                                      min_cells_per_group = 50,
                                      min_counts = 10,
                                      robust = TRUE,
                                      seed = 42,
                                      offset_vector_name = NULL,
                                      offset_vector_as = c("od","ratio")) {
  offset_vector_as <- match.arg(offset_vector_as)
  results <- list()
  
  sample_names <- names(seurat_list)
  if (is.null(sample_names)) sample_names <- paste0("sample_", seq_along(seurat_list))
  
  # Pre-extract feature vector if requested (same across samples)
  feature_vec <- NULL
  if (!is.null(offset_vector_name)) {
    feature_vec <- extract_feature_vector(seurat_list, offset_vector_name)
  }
  
  for (assay in assay_names) {
    .log("\n--- Calculating within-sample BCV for %s assay ---", assay)
    assay_results <- list()
    
    for (i in seq_along(seurat_list)) {
      sample_name <- sample_names[i]
      obj <- seurat_list[[i]]
      if (!(assay %in% names(obj@assays))) {
        warning(sprintf("Assay %s not found in %s, skipping", assay, sample_name))
        next
      }
      .log("Processing %s (%d cells)", sample_name, ncol(obj))
      tryCatch({
        bcv_result <- calculate_within_sample_bcv(
          obj, assay_name = assay,
          n_groups = n_groups,
          cells_per_group = cells_per_group,
          min_cells_per_group = min_cells_per_group,
          min_counts = min_counts,
          robust = robust,
          seed = seed + i,
          feature_vector = feature_vec,   # NULL if not requested
          vector_as = offset_vector_as
        )
        bcv_result$sample_name <- sample_name
        assay_results[[sample_name]] <- bcv_result
      }, error = function(e) {
        warning(sprintf("Error processing %s for %s: %s", sample_name, assay, e$message))
      })
    }
    if (length(assay_results) > 0) results[[assay]] <- assay_results
  }
  results
}

#' Summarize Within-Sample BCV Results
#'
#' Collapses a list of detailed within-sample BCV results (from
#' \code{\link{analyze_within_sample_bcv}}) into a summary data frame.
#'
#' @param within_bcv_results The nested list returned by
#'   \code{\link{analyze_within_sample_bcv}}.
#'
#' @return A `data.frame` summarizing the BCV metrics for each sample and
#'   assay. Key columns include:
#'   \itemize{
#'     \item `sample`: The sample identifier.
#'     \item `Assay`: The assay name (e.g., "RNA", "RNA_corr").
#'     \item `N_Cells`, `N_Groups`, `N_Genes`: Counts of cells, groups, and genes.
#'     \item `Common_BCV`, `Median_BCV`, `Mean_BCV`, `SD_BCV`: Summary statistics
#'       for the biological coefficient of variation.
#'     \item `Reduction_Percent`: If both "RNA" and "RNA_corr" assays are present,
#'       this column shows the percentage reduction in median BCV from RNA to
#'       RNA_corr.
#'   }
#'
#' @export
summarize_within_sample_bcv <- function(within_bcv_results) {
  summary_list <- list()
  for (assay in names(within_bcv_results)) {
    assay_data <- within_bcv_results[[assay]]
    sample_stats <- lapply(names(assay_data), function(sample) {
      res <- assay_data[[sample]]
      data.frame(
        sample = sample,
        Assay = assay,
        N_Cells = res$n_cells_total,
        N_Groups = res$n_groups,
        N_Genes = res$n_genes,
        Common_BCV = res$bcv_common,
        Median_BCV = median(res$bcv_tagwise, na.rm = TRUE),
        Mean_BCV = mean(res$bcv_tagwise, na.rm = TRUE),
        SD_BCV = sd(res$bcv_tagwise, na.rm = TRUE)
      )
    })
    summary_list[[assay]] <- do.call(rbind, sample_stats)
  }
  summary_df <- do.call(rbind, summary_list)
  rownames(summary_df) <- NULL
  
  if (all(c("RNA", "RNA_corr") %in% names(within_bcv_results))) {
    rna_df <- summary_list[["RNA"]]
    corr_df <- summary_list[["RNA_corr"]]
    merged <- merge(rna_df[, c("sample", "Median_BCV")], 
                    corr_df[, c("sample", "Median_BCV")],
                    by = "sample", suffixes = c("_RNA", "_RNA_corr"))
    merged$Reduction_Percent <- round(
      (merged$Median_BCV_RNA - merged$Median_BCV_RNA_corr) / merged$Median_BCV_RNA * 100, 2
    )
    summary_df <- merge(summary_df, merged[, c("sample", "Reduction_Percent")],
                        by = "sample", all.x = TRUE)
    summary_df$Reduction_Percent[is.na(summary_df$Reduction_Percent)] <- 0
  }
  summary_df
}



#' @name plot_within_sample_bcv
#' @title Plot within-sample BCV for a specific sample (or show a summary)
#'
#' @description
#' Visualizes within-sample BCV (tagwise) vs average log2 CPM for a selected
#' sample from \code{\link{analyze_within_sample_bcv}}/\code{\link{calculate_within_sample_bcv}}
#' results. If \code{sample_name} is \code{NULL}, returns the overall summary
#' from \code{\link{plot_within_sample_summary}}.
#'
#' @param within_bcv_results List of BCV results (per assay) as returned by
#'   \code{\link{analyze_within_sample_bcv}}.
#' @param sample_name Optional character; the sample to plot. If \code{NULL},
#'   a summary plot is returned instead.
#'
#' @return A \code{ggplot} object (or patchwork of two plots if both assays are
#'   present) for the requested sample, or the summary plot if \code{NULL}.
#'
#' @seealso \code{\link{analyze_within_sample_bcv}},
#'   \code{\link{plot_within_sample_summary}}
#' @export


plot_within_sample_bcv <- function(within_bcv_results, sample_name = NULL) {
  if (is.null(sample_name)) return(plot_within_sample_summary(within_bcv_results))
  plots <- list()
  for (assay in names(within_bcv_results)) {
    if (sample_name %in% names(within_bcv_results[[assay]])) {
      res <- within_bcv_results[[assay]][[sample_name]]
      df <- data.frame(ave_log_cpm = res$ave_log_cpm, bcv = res$bcv_tagwise)
      p <- ggplot(df, aes(x = ave_log_cpm, y = bcv)) +
        geom_point(alpha = 0.5, size = 0.8, color = ifelse(assay == "RNA", "darkblue", "darkgreen")) +
        geom_line(data = data.frame(ave_log_cpm = res$ave_log_cpm, bcv = res$bcv_trend),
                  aes(y = bcv), color = "orange", size = 1.2) +
        geom_hline(yintercept = res$bcv_common, linetype = "dashed", color = "red", size = 1) +
        labs(x = "Average log2 CPM", y = "BCV (within-sample)",
             title = sprintf("%s - %s (within-sample)", sample_name, assay),
             subtitle = sprintf("Common BCV: %.3f, Median: %.3f", res$bcv_common,
                                median(res$bcv_tagwise, na.rm = TRUE))) +
        theme_classic()
      plots[[assay]] <- p
    }
  }
  if (length(plots) == 2) return(plots$RNA + plots$RNA_corr)
  if (length(plots) == 1) return(plots[[1]])
  NULL
}


#' Plot a Summary of Within-Sample BCV Results
#'
#' Creates summary visualizations for within-sample BCV results across multiple
#' samples and assays.
#'
#' @param within_bcv_results A nested list of BCV results as returned by
#'   \code{\link{analyze_within_sample_bcv}}.
#'
#' @return A `ggplot` object, which is a composite of two plots if both "RNA" and
#' "RNA_corr" assays are present:
#' \itemize{
#'   \item A bar plot comparing the median within-sample BCV for each assay
#'     across all samples.
#'   \item A scatter plot directly comparing the median BCV of the "RNA" assay
#'     against the "RNA_corr" assay for each sample.
#' }
#' If only one assay is present, only the bar plot is returned. Returns an empty
#' plot with a message if no valid results are provided.
#'
#' @details This function provides a high-level overview of the impact of
#' overdispersion correction on within-sample variability. The scatter plot is
#' particularly useful for assessing the consistency of BCV reduction across
#' different samples.
#'
#' @seealso \code{\link{analyze_within_sample_bcv}}, \code{\link{plot_within_sample_bcv}}
#' @export


plot_within_sample_summary <- function(within_bcv_results) {
  plot_data <- list()
  for (assay in names(within_bcv_results)) {
    for (sample in names(within_bcv_results[[assay]])) {
      res <- within_bcv_results[[assay]][[sample]]
      if (is.null(res)) next
      plot_data <- append(plot_data, list(data.frame(
        sample = sample, Assay = assay,
        Common_BCV = res$bcv_common,
        Median_BCV = median(res$bcv_tagwise, na.rm = TRUE),
        Mean_BCV = mean(res$bcv_tagwise, na.rm = TRUE)
      )))
    }
  }
  if (!length(plot_data)) {
    return(
      ggplot() +
        annotate("text", x = 0.5, y = 0.5, label = "No within-sample BCV results available", size = 5) +
        theme_void()
    )
  }
  df <- do.call(rbind, plot_data)
  
  p1 <- ggplot(df, aes(x = sample, y = Median_BCV, fill = Assay)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(y = "Median BCV (within-sample)", title = "Within-sample BCV across samples") +
    theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  if (all(c("RNA", "RNA_corr") %in% unique(df$Assay))) {
    df_wide <- df %>% select(sample, Assay, Median_BCV) %>%
      tidyr::pivot_wider(names_from = Assay, values_from = Median_BCV)
    p2 <- ggplot(df_wide, aes(x = RNA, y = RNA_corr)) +
      geom_point(size = 3) + geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
      geom_text(aes(label = sample), hjust = -0.1, vjust = -0.1, size = 3) +
      labs(x = "Median BCV - RNA", y = "Median BCV - RNA_corr", title = "Within-sample BCV comparison") +
      theme_classic() + coord_equal()
    return(p1 / p2)
  }
  p1
}


#' Extract Common Features and Create a Pseudobulk Matrix Across Samples
#'
#' This function identifies features (e.g., genes) that are common across all 
#' \code{Seurat} objects in a list, extracts their counts from a specified assay, 
#' and aggregates counts across cells within each object to create a pseudobulk 
#' matrix.
#'
#' @param seurat_list A list of \code{Seurat} objects containing the assay of interest.
#' @param assay_name Character scalar specifying the assay from which to extract counts 
#'   (default: \code{"RNA"}).
#'
#' @return A sparse matrix (\code{dgCMatrix}) of pseudobulk counts with features 
#'   as rows and samples as columns. Column names correspond to \code{names(seurat_list)}, 
#'   or are automatically generated if missing.
#'
#' @details
#' The function:
#' \enumerate{
#'   \item Finds the intersection of feature names across all \code{Seurat} objects.
#'   \item Extracts the raw counts for these common features from the specified assay.
#'   \item Aggregates counts across all cells within each object to produce one pseudobulk 
#'         column per sample.
#' }
#'
#' @examples
#' \dontrun{
#' pb <- extract_and_pseudobulk(my_seurat_list, assay_name = "RNA")
#' dim(pb) # features x samples
#' }
#'
#' @export


extract_and_pseudobulk <- function(seurat_list, assay_name = "RNA") {
  .log("Processing %d objects for %s assay", length(seurat_list), assay_name)
  all_features <- lapply(seurat_list, function(obj) rownames(obj[[assay_name]]))
  common_features <- Reduce(intersect, all_features)
  .log("Found %d common features", length(common_features))
  
  nm <- names(seurat_list); if (is.null(nm) || any(is.na(nm) | nm == "")) nm <- paste0("sample_", seq_along(seurat_list))
  
  pb_cols <- vector("list", length(seurat_list))
  for (i in seq_along(seurat_list)) {
    .log("  %s: %d cells", nm[i], ncol(seurat_list[[i]]))
    counts <- GetAssayData(seurat_list[[i]], assay = assay_name, layer = "counts")
    counts <- counts[common_features, , drop = FALSE]
    pb_cols[[i]] <- Matrix::rowSums(counts)
    rm(counts); gc(FALSE)
  }
  pb <- do.call(cbind, pb_cols)
  rownames(pb) <- common_features
  colnames(pb) <- nm
  as(pb, "dgCMatrix")
}


#' Diagnose presence and effect of a corrected assay in a Seurat object
#'
#' Prints a report comparing totals and gene-level ratios between the raw "RNA"
#' and corrected "RNA_corr" assays, and shows basic overdispersion metadata if
#' available.
#'
#' @param objs List of \code{Seurat} objects.
#' @param genes_to_check Optional character vector of genes to inspect; if
#'   \code{NULL}, representative genes are selected by expression quantiles.
#' @param sample_idx Integer; which object to inspect (default 1).
#'
#' @return Invisibly returns a list with \code{gene_ratios} and \code{cell_ratios}.
#'
#' @seealso \code{\link{process_and_create_seurat_corrected_improved}}
#' @export



diagnose_seurat_correction <- function(objs, genes_to_check = NULL, sample_idx = 1) {
  obj <- objs[[sample_idx]]
  sample_name <- names(objs)[sample_idx]; if (is.null(sample_name)) sample_name <- paste0("sample_", sample_idx)
  cat("\n", strrep("=", 60), "\n")
  cat(sprintf("DIAGNOSTIC REPORT FOR: %s\n", sample_name))
  cat(strrep("=", 60), "\n")
  cat("\nAvailable assays:", paste(names(obj@assays), collapse = ", "), "\n")
  if (!all(c("RNA","RNA_corr") %in% names(obj@assays))) { cat("ERROR: Both RNA and RNA_corr assays are required\n"); return(invisible(NULL)) }
  raw_counts  <- GetAssayData(obj, assay = "RNA",      layer = "counts")
  corr_counts <- GetAssayData(obj, assay = "RNA_corr", layer = "counts")
  feature_meta <- obj@assays[["RNA_corr"]]@misc[["feature_meta"]]
  if (!is.null(feature_meta)) {
    cat("\nFeature metadata columns in RNA_corr:", paste(colnames(feature_meta), collapse = ", "), "\n")
    if ("OverDisp" %in% names(feature_meta)) {
      od_stats <- summary(feature_meta$OverDisp)
      cat("\nOverdispersion statistics:\n"); print(od_stats)
      cat(sprintf("  Transcripts with OD > 1: %d (%.1f%%)\n", sum(feature_meta$OverDisp > 1, na.rm = TRUE),
                  100 * mean(feature_meta$OverDisp > 1, na.rm = TRUE)))
    }
  }
  common_genes <- intersect(rownames(raw_counts), rownames(corr_counts))
  raw_counts  <- raw_counts[common_genes, ]
  corr_counts <- corr_counts[common_genes, ]
  total_raw  <- Matrix::colSums(raw_counts)
  total_corr <- Matrix::colSums(corr_counts)
  ratio_per_cell <- total_corr / total_raw
  cat("\n--- Per-cell totals ---\n")
  cat(sprintf("Median total counts (raw): %.0f\n", median(total_raw)))
  cat(sprintf("Median total counts (corrected): %.0f\n", median(total_corr)))
  cat(sprintf("Median ratio (corr/raw): %.3f\n", median(ratio_per_cell)))
  if (is.null(genes_to_check)) {
    gene_means <- Matrix::rowMeans(raw_counts)
    q <- quantile(gene_means[gene_means > 0], c(0.1, 0.5, 0.9))
    genes_to_check <- c(
      names(gene_means)[which.min(abs(gene_means - q[1]))],
      names(gene_means)[which.min(abs(gene_means - q[2]))],
      names(gene_means)[which.min(abs(gene_means - q[3]))]
    )
  }
  cat("\n--- Gene-level inspection ---\n")
  for (gene in genes_to_check) {
    if (!(gene %in% common_genes)) next
    raw_g  <- as.numeric(raw_counts[gene, ])
    corr_g <- as.numeric(corr_counts[gene, ])
    nonzero <- raw_g > 0 & corr_g > 0
    if (sum(nonzero) > 0) {
      ratios <- corr_g[nonzero] / raw_g[nonzero]
      cat(sprintf("\nGene: %s\n", gene))
      cat(sprintf("  Cells with expression: %d/%d\n", sum(nonzero), length(nonzero)))
      cat(sprintf("  Mean raw count: %.2f\n", mean(raw_g[nonzero])))
      cat(sprintf("  Mean corrected count: %.2f\n", mean(corr_g[nonzero])))
      cat(sprintf("  Median ratio (corr/raw): %.3f\n", median(ratios)))
      if (!is.null(feature_meta) && "OverDisp" %in% names(feature_meta)) {
        if (gene %in% rownames(feature_meta)) {
          cat(sprintf("  Overdispersion: %.3f\n", feature_meta[gene, "OverDisp"]))
        }
      }
    }
  }
  cat("\n--- Overall gene statistics ---\n")
  gene_sums_raw  <- Matrix::rowSums(raw_counts)
  gene_sums_corr <- Matrix::rowSums(corr_counts)
  gene_ratios <- gene_sums_corr / gene_sums_raw
  gene_ratios <- gene_ratios[is.finite(gene_ratios) & gene_ratios > 0]
  cat(sprintf("Median gene ratio (corr/raw): %.3f\n", median(gene_ratios)))
  cat(sprintf("Mean gene ratio: %.3f\n", mean(gene_ratios)))
  cat(sprintf("Transcripts with ratio > 1: %d (%.1f%%)\n", sum(gene_ratios > 1), 100 * mean(gene_ratios > 1)))
  cat("\n", strrep("=", 60), "\n\n")
  invisible(list(
    gene_ratios = gene_ratios,
    cell_ratios = ratio_per_cell
  ))
}



#' Compute BCV directly from pseudobulk counts
#'
#' Estimates BCV from a gene-by-sample pseudobulk matrix using edgeR. This is
#' the simplest between-sample BCV baseline (no additional offsets).
#'
#' @param pb_counts Gene-by-sample count matrix (\code{matrix} or \code{dgCMatrix}).
#' @param min_counts Integer; minimum total counts per gene to retain (default 10).
#' @param robust Logical; robust estimation in \code{edgeR::estimateDisp} (default TRUE).
#' @param filter_method Either \code{"rowSums"} or \code{"filterByExpr"} (default "rowSums").
#'
#' @return A list with fields \code{genes}, \code{samples}, \code{ave_log_cpm},
#'   \code{bcv_tagwise}, \code{bcv_common}, \code{bcv_trend}, \code{n_genes},
#'   and \code{n_samples}.
#'
#' @seealso \code{\link{calculate_bcv}}, \code{\link{plot_bcv}},
#'   \code{\link{plot_bcv_comparison}}, \code{\link{create_summary_table}}
#' @export


calculate_bcv_direct <- function(pb_counts, min_counts = 10, robust = TRUE,
                                 filter_method = c("rowSums","filterByExpr")) {
  filter_method <- match.arg(filter_method)
  if (ncol(pb_counts) < 2) stop("Need at least 2 samples to calculate BCV")
  group <- rep(1, ncol(pb_counts))
  keep <- if (filter_method == "rowSums") rowSums(pb_counts) >= min_counts
  else edgeR::filterByExpr(pb_counts, group = group, min.count = min_counts)
  if (!any(keep)) stop("[calculate_bcv_direct] No genes passed filtering.")
  y_counts <- pb_counts[keep, , drop = FALSE]
  .log("Kept %d/%d genes after filtering", sum(keep), length(keep))
  y <- edgeR::DGEList(counts = y_counts)
  y <- edgeR::calcNormFactors(y, method = "TMM")
  design <- stats::model.matrix(~ 1, data = data.frame(sample = colnames(y_counts)))
  y <- edgeR::estimateDisp(y, design = design, robust = robust)
  list(
    genes = rownames(y$counts),
    samples = colnames(y$counts),
    ave_log_cpm = y$AveLogCPM,
    bcv_tagwise = sqrt(y$tagwise.dispersion),
    bcv_common = sqrt(y$common.dispersion),
    bcv_trend = sqrt(y$trended.dispersion),
    n_genes = nrow(y$counts),
    n_samples = ncol(y$counts)
  )
}



#' Calculate BCV with Optional Correction Offsets
#'
#' Computes the Biological Coefficient of Variation (BCV) from a raw pseudobulk
#' count matrix, with the option to apply gene-specific offsets to model
#' technical variability (e.g., from overdispersion).
#'
#' @param pb_raw A raw pseudobulk count matrix (genes-by-samples).
#' @param pb_corr An optional corrected pseudobulk count matrix. Required if
#'   `offset_method = "ratio"`.
#' @param overdispersion Deprecated. Use `feature_vector` instead. A named numeric
#'   vector of per-gene overdispersion factors.
#' @param feature_vector A named numeric vector of per-gene values to be used as
#'   an offset (e.g., overdispersion factors or scaling ratios).
#' @param vector_as How to interpret `feature_vector`: `"od"` (overdispersion) or
#'   `"ratio"` (scaling factor).
#' @param min_counts Minimum total count for a gene to be included in the analysis.
#' @param robust Logical; whether to use robust dispersion estimation in `edgeR`.
#' @param filter_method Method for gene filtering: `"rowSums"` or `"filterByExpr"`.
#' @param offset_method The method for applying corrections:
#'   \describe{
#'     \item{`"ratio"`}{Calculates offsets from the log-ratio of corrected to raw counts (`log(pb_corr / pb_raw)`).}
#'     \item{`"vector"`}{Calculates offsets from the provided `feature_vector`.}
#'     \item{`"overdispersion"`}{Legacy alias for `"vector"` when using an overdispersion vector.}
#'   }
#'
#' @details This function serves as the core of between-sample BCV analysis when
#' accounting for technical noise. By providing a `feature_vector` (like
#' scaling_factor), the `edgeR` model can distinguish biological
#' variability from known technical variance, ideally resulting in a lower and
#' more accurate BCV estimate.
#'
#' @return A list containing standard `edgeR` dispersion results, including
#' `bcv_tagwise`, `bcv_common`, `ave_log_cpm`, and other metrics.
#'
#' @examples
#' \dontrun{
#' # Assuming pb_raw is a raw counts matrix and od_vector contains per-gene OD
#' bcv_results_corrected <- calculate_bcv(
#'   pb_raw,
#'   feature_vector = od_vector,
#'   offset_method = "vector",
#'   vector_as = "od"
#' )
#' }
#' @export


calculate_bcv <- function(pb_raw,
                          pb_corr = NULL,
                          overdispersion = NULL,   # legacy path; kept for compatibility
                          feature_vector = NULL,
                          vector_as = c("od","ratio"),
                          min_counts = 10,
                          robust = TRUE,
                          filter_method = c("rowSums","filterByExpr"),
                          offset_method = c("ratio", "vector", "overdispersion")) {
  filter_method <- match.arg(filter_method)
  offset_method <- match.arg(offset_method)
  vector_as <- match.arg(vector_as)
  if (ncol(pb_raw) < 2) stop("Need at least 2 samples to calculate BCV")
  
  group <- rep(1, ncol(pb_raw))
  keep <- if (filter_method == "rowSums") rowSums(pb_raw) >= min_counts
  else edgeR::filterByExpr(pb_raw, group = group, min.count = min_counts)
  if (!any(keep)) stop("[calculate_bcv] No genes passed filtering.")
  y_counts <- pb_raw[keep, , drop = FALSE]
  .log("Kept %d/%d genes after filtering", sum(keep), length(keep))
  
  y <- edgeR::DGEList(counts = y_counts)
  y <- edgeR::calcNormFactors(y, method = "TMM")
  design <- stats::model.matrix(~ 1, data = data.frame(sample = colnames(y_counts)))
  
  raw_mat <- as.matrix(y_counts)
  base_off <- matrix(
    log(y$samples$lib.size * y$samples$norm.factors),
    nrow = nrow(raw_mat), ncol = ncol(raw_mat), byrow = TRUE
  )
  
  if (offset_method == "ratio") {
    if (is.null(pb_corr)) stop("offset_method='ratio' requires pb_corr")
    pb_corr_keep <- pb_corr[rownames(y_counts), colnames(y_counts), drop = FALSE]
    corr_mat <- as.matrix(pb_corr_keep)
    corr_mat[corr_mat < 1e-12] <- 1e-12
    raw_mat[raw_mat < 1e-12] <- 1e-12
    K_eff <- corr_mat / raw_mat                # corrected/raw
    corr_off <- log(K_eff); corr_off[!is.finite(corr_off)] <- 0
    y$offset <- base_off + corr_off
    .log("Applied ratio-based offsets (median K_eff = %.3f)", median(K_eff, na.rm = TRUE))
    storage.mode(y$offset) <- "double"
    
  } else if (offset_method == "vector") {
    if (is.null(feature_vector)) stop("offset_method='vector' requires feature_vector")
    fv <- feature_vector[rownames(y_counts)]
    if (any(is.na(fv))) {
      warning("Some genes missing feature_vector values, setting to 1")
      fv[is.na(fv)] <- 1
    }
    fv <- pmax(fv, 1e-12)
    add_off <- matrix(log(fv), nrow = length(fv), ncol = ncol(raw_mat), byrow = FALSE)
    y$offset <- base_off + add_off
    .log("Applied vector-based offsets using feature vector interpreted as '%s'", vector_as)
    storage.mode(y$offset) <- "double"
    
  } else if (offset_method == "overdispersion") {  # legacy alias
    if (is.null(overdispersion)) stop("offset_method='overdispersion' requires overdispersion vector")
    od_filtered <- overdispersion[rownames(y_counts)]
    if (any(is.na(od_filtered))) {
      warning("Some genes missing overdispersion values, setting to 1")
      od_filtered[is.na(od_filtered)] <- 1
    }
    od_filtered <- pmax(od_filtered, 1)
    od_off <- matrix(log(od_filtered), nrow = length(od_filtered), ncol = ncol(raw_mat), byrow = FALSE)
    y$offset <- base_off + od_off
    .log("Applied overdispersion-based offsets (median OD = %.3f)", median(od_filtered, na.rm = TRUE))
    storage.mode(y$offset) <- "double"
  }
  
  y <- edgeR::estimateDisp(y, design = design, robust = robust)
  list(
    genes = rownames(y$counts),
    samples = colnames(y$counts),
    ave_log_cpm = y$AveLogCPM,
    bcv_tagwise = sqrt(y$tagwise.dispersion),
    bcv_common = sqrt(y$common.dispersion),
    bcv_trend = sqrt(y$trended.dispersion),
    n_genes = nrow(y$counts),
    n_samples = ncol(y$counts)
  )
}


#' Plot Biological Coefficient of Variation (BCV)
#'
#' Creates a scatter plot of tagwise BCV versus average log2 counts per million (CPM),
#' with overlaid trended BCV and the common BCV line.
#'
#' @param bcv_results A list returned by \code{\link{calculate_bcv_direct}} or
#'   \code{\link{calculate_within_sample_bcv}}, containing \code{ave_log_cpm},
#'   \code{bcv_tagwise}, \code{bcv_common}, and \code{bcv_trend}.
#' @param title Character; the plot title (default: \code{"BCV Plot"}).
#' @param color Character; color of the scatter points (default: \code{"darkblue"}).
#'
#' @return A \pkg{ggplot2} object showing:
#' \itemize{
#'   \item Points for tagwise BCV vs. average log2 CPM.
#'   \item An orange trend line for trended BCV.
#'   \item A dashed red horizontal line for the common BCV.
#' }
#'
#' @details
#' This plot is useful for visualizing the mean–variance relationship in
#' RNA-seq data and checking the fit of dispersion models.
#'
#' @examples
#' \dontrun{
#' pb <- matrix(rpois(1000, lambda = 10), nrow = 100, ncol = 10)
#' res <- calculate_bcv_direct(pb)
#' plot_bcv(res, title = "Example BCV Plot")
#' }
#'
#' @export


plot_bcv <- function(bcv_results, title = "BCV Plot", color = "darkblue") {
  df <- data.frame(ave_log_cpm = bcv_results$ave_log_cpm, bcv = bcv_results$bcv_tagwise)
  common_bcv <- bcv_results$bcv_common
  trend_df <- data.frame(ave_log_cpm = bcv_results$ave_log_cpm, bcv_trend = bcv_results$bcv_trend)
  ggplot(df, aes(x = ave_log_cpm, y = bcv)) +
    geom_point(alpha = 0.5, color = color, size = 0.8) +
    geom_line(data = trend_df, aes(y = bcv_trend), color = "orange", size = 1.2) +
    geom_hline(yintercept = common_bcv, linetype = "dashed", color = "red", size = 1) +
    labs(x = "Average log2 CPM", y = "Biological Coefficient of Variation (BCV)", title = title) +
    theme_classic() +
    coord_cartesian(ylim = c(0, max(bcv_results$bcv_tagwise) * 1.05))
}

#' Compare BCV Between Raw and Corrected Assays
#'
#' Creates a scatter plot comparing tagwise Biological Coefficient of Variation (BCV)
#' between raw (\code{RNA}) and corrected (\code{RNA_corr}) assays for common genes.
#'
#' @param bcv_rna A BCV results list (from \code{\link{calculate_bcv_direct}} or
#'   \code{\link{calculate_within_sample_bcv}}) for the raw \code{RNA} assay.
#' @param bcv_corr A BCV results list for the corrected \code{RNA_corr} assay.
#'
#' @return A \pkg{ggplot2} object showing:
#' \itemize{
#'   \item Points representing per-gene BCV values in \code{RNA} (x-axis) vs.
#'         \code{RNA_corr} (y-axis), colored by average log2 CPM.
#'   \item A dashed red identity line for reference.
#'   \item A title reporting the Pearson correlation and median BCV reduction percentage.
#' }
#'
#' @details
#' The plot highlights how overdispersion correction affects BCV across genes.
#' The color gradient corresponds to average log2 CPM, providing insight into
#' whether effects differ by gene abundance.
#'
#' @examples
#' \dontrun{
#' # Assuming bcv_rna and bcv_corr are BCV result lists for raw and corrected assays
#' plot_bcv_comparison(bcv_rna, bcv_corr)
#' }
#'
#' @export


plot_bcv_comparison <- function(bcv_rna, bcv_corr) {
  common <- intersect(bcv_rna$genes, bcv_corr$genes)
  if (!length(common)) return(ggplot() + labs(title = "No common genes found") + theme_classic())
  idx_r <- match(common, bcv_rna$genes)
  idx_c <- match(common, bcv_corr$genes)
  df <- data.frame(
    bcv_rna  = bcv_rna$bcv_tagwise[idx_r],
    bcv_corr = bcv_corr$bcv_tagwise[idx_c],
    ave_log_cpm = bcv_rna$ave_log_cpm[idx_r]
  )
  cor_val <- cor(df$bcv_rna, df$bcv_corr, use = "complete.obs")
  med_red <- (median(df$bcv_rna, na.rm = TRUE) - median(df$bcv_corr, na.rm = TRUE)) /
    median(df$bcv_rna, na.rm = TRUE) * 100
  ggplot(df, aes(x = bcv_rna, y = bcv_corr)) +
    geom_point(aes(color = ave_log_cpm), alpha = 0.5, size = 0.8) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
    scale_color_gradient2(low = "blue", mid = "gray", high = "red",
                          midpoint = median(df$ave_log_cpm), name = "Ave log CPM") +
    labs(x = "BCV (RNA)", y = "BCV (RNA_corr)",
         title = sprintf("BCV Comparison (r = %.3f, reduction = %.1f%%)", cor_val, med_red)) +
    theme_classic() + coord_equal()
}


#' Create a between/within BCV comparison plot
#'
#' Compares the distribution of BCV values between samples and within samples
#' (across assays), using violin/boxplot summaries and bar summaries of medians.
#'
#' @param between_results A named list of between-sample BCV results, typically
#'   containing elements \code{$RNA} and (optionally) \code{$RNA_corr}, each as
#'   returned by \code{\link{calculate_bcv_direct}} or \code{\link{calculate_bcv}}.
#' @param within_results The within-sample BCV results list returned by
#'   \code{\link{analyze_within_sample_bcv}}.
#'
#' @return A \code{ggplot}/patchwork object showing distributions and medians.
#'
#' @seealso \code{\link{calculate_bcv_direct}}, \code{\link{calculate_bcv}},
#'   \code{\link{analyze_within_sample_bcv}}
#' @export


plot_between_vs_within_bcv <- function(between_results, within_results) {
  between_data <- lapply(names(between_results), function(assay)
    data.frame(Type = "Between-sample", Assay = assay, BCV = between_results[[assay]]$bcv_tagwise))
  all_between <- do.call(rbind, between_data)
  
  # Build within only if present
  within_data <- list()
  for (assay in names(within_results)) {
    for (sample in names(within_results[[assay]])) {
      res <- within_results[[assay]][[sample]]
      if (is.null(res)) next
      within_data[[paste(assay, sample, sep = "_")]] <- data.frame(
        Type = "Within-sample", Assay = assay, BCV = res$bcv_tagwise, sample = sample
      )
    }
  }
  all_within <- if (length(within_data)) do.call(rbind, within_data) else NULL
  if (!is.null(all_within)) all_within$sample <- NULL
  
  all_data <- if (is.null(all_within)) all_between else rbind(all_between, all_within)
  
  p1 <- ggplot(all_data, aes(x = interaction(Type, Assay), y = BCV, fill = Assay)) +
    geom_violin(trim = FALSE, alpha = 0.7) +
    geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.5) +
    labs(x = "", y = "Biological Coefficient of Variation", title = "Between-sample vs Within-sample BCV") +
    theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_manual(values = c("RNA" = "darkblue", "RNA_corr" = "darkgreen"))
  
  summary_stats <- all_data %>% dplyr::group_by(Type, Assay) %>%
    dplyr::summarise(Median_BCV = median(BCV, na.rm = TRUE),
                     Mean_BCV = mean(BCV, na.rm = TRUE), .groups = "drop")
  
  if (!nrow(summary_stats)) {
    return(p1)
  }
  
  p2 <- ggplot(summary_stats, aes(x = Type, y = Median_BCV, fill = Assay)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(y = "Median BCV", title = "Median BCV Comparison") +
    theme_classic() +
    scale_fill_manual(values = c("RNA" = "darkblue", "RNA_corr" = "darkgreen"))
  
  p1 / p2
}


#' Create a Summary Table for Between-Sample BCV Results
#'
#' Generates a concise data frame summarizing the key metrics from a list of
#' between-sample BCV analysis results.
#'
#' @param results A named list where each element is a BCV result object as
#'   returned by `calculate_bcv` or `calculate_bcv_direct`. The names of the
#'   list elements (e.g., "RNA", "RNA_corr") are used to identify the assay.
#'
#' @return A `data.frame` with one row per assay, containing columns for:
#'   \itemize{
#'     \item `Assay`: The name of the assay.
#'     \item `N_Genes`, `N_samples`: Number of genes and samples in the analysis.
#'     \item `Common_BCV`, `Median_BCV`, `Mean_BCV`, `Min_BCV`, `Max_BCV`:
#'       Summary statistics for the biological coefficient of variation.
#'     \item `Reduction_Percent`: If both "RNA" and "RNA_corr" are present, this
#'       column shows the percentage reduction in median BCV for the corrected assay.
#'   }
#'
#' @seealso \code{\link{calculate_bcv}}, \code{\link{print_summary}}
#' @export
create_summary_table <- function(results) {
  lst <- lapply(names(results), function(an) {
    b <- results[[an]]
    data.frame(
      Assay = an, 
      N_Genes = b$n_genes, 
      N_samples = b$n_samples,
      Common_BCV = round(b$bcv_common, 4),
      Median_BCV = round(median(b$bcv_tagwise, na.rm = TRUE), 4),
      Mean_BCV = round(mean(b$bcv_tagwise, na.rm = TRUE), 4),
      Min_BCV = round(min(b$bcv_tagwise, na.rm = TRUE), 4),
      Max_BCV = round(max(b$bcv_tagwise, na.rm = TRUE), 4)
    )
  })
  out <- do.call(rbind, lst); rownames(out) <- NULL
  if (all(c("RNA","RNA_corr") %in% names(results))) {
    mr <- out$Median_BCV[out$Assay=="RNA"]
    mc <- out$Median_BCV[out$Assay=="RNA_corr"]
    out$Reduction_Percent <- ifelse(out$Assay=="RNA_corr", round((mr - mc)/mr*100, 2), 0)
  }
  out
}


#' Print Summary of BCV Analysis Results
#'
#' Displays a formatted summary of Biological Coefficient of Variation (BCV)
#' analysis results for either between-sample or within-sample comparisons.
#'
#' @param results A list of BCV results.
#'   \itemize{
#'     \item For \code{type = "between"}: each element is a BCV result object
#'       (as from \code{\link{calculate_bcv_direct}}).
#'     \item For \code{type = "within"}: each element is itself a list of
#'       BCV results per sample (as from \code{\link{calculate_within_sample_bcv}}).
#'   }
#' @param type Character string indicating whether the results are from a
#'   \code{"between"}-sample analysis or a \code{"within"}-sample analysis.
#'   Defaults to \code{"between"}.
#'
#' @details
#' The summary includes the number of genes, samples, common BCV, median
#' tagwise BCV, and BCV ranges. If both \code{RNA} and \code{RNA_corr} are
#' present, it also reports the percent reduction in BCV for \code{RNA_corr}
#' relative to \code{RNA}.
#'
#' @return Invisibly returns \code{NULL}. Output is printed to the console.
#'
#' @examples
#' \dontrun{
#' print_summary(between_results, type = "between")
#' print_summary(within_results, type = "within")
#' }
#'
#' @seealso
#' \code{\link{calculate_bcv_direct}},
#' \code{\link{calculate_within_sample_bcv}},
#' \code{\link{create_summary_table}}
#'
#' @export


print_summary <- function(results, type = "between") {
  cat("\n", strrep("=", 60), "\n"); cat(sprintf("BCV ANALYSIS SUMMARY (%s-sample)\n", type))
  cat(strrep("=", 60), "\n\n")
  for (nm in names(results)) {
    b <- results[[nm]]
    cat(sprintf(">>> %s:\n", nm))
    if (type == "between") {
      cat(sprintf("    samples: %d\n", b$n_samples))
      cat(sprintf("    Transcripts analyzed: %d\n", b$n_genes))
      cat(sprintf("    Common BCV: %.4f\n", b$bcv_common))
      cat(sprintf("    Median tagwise BCV: %.4f\n", median(b$bcv_tagwise, na.rm = TRUE)))
      cat(sprintf("    BCV range: [%.4f, %.4f]\n\n",
                  min(b$bcv_tagwise, na.rm = TRUE), max(b$bcv_tagwise, na.rm = TRUE)))
    } else {
      sample_medians <- sapply(b, function(x) median(x$bcv_tagwise, na.rm = TRUE))
      cat(sprintf("    Number of samples analyzed: %d\n", length(b)))
      cat(sprintf("    Median BCV across samples: %.4f\n", median(sample_medians)))
      cat(sprintf("    BCV range across samples: [%.4f, %.4f]\n\n", min(sample_medians), max(sample_medians)))
    }
  }
  if (all(c("RNA","RNA_corr") %in% names(results))) {
    if (type == "between") {
      mr <- median(results$RNA$bcv_tagwise, na.rm = TRUE)
      mc <- median(results$RNA_corr$bcv_tagwise, na.rm = TRUE)
      red <- (mr - mc)/mr*100
      cat(strrep("-", 60), "\n"); cat("COMPARISON:\n")
      cat(sprintf("  Median BCV reduction: %.2f%%\n", red))
      cat(sprintf("  Common BCV reduction: %.2f%%\n",
                  (results$RNA$bcv_common - results$RNA_corr$bcv_common)/
                    results$RNA$bcv_common*100))
      if (red > 0) cat("  -> RNA_corr shows REDUCED biological variation \n")
      else cat("  -> RNA_corr shows INCREASED biological variation \n")
    } else {
      rna_medians  <- sapply(results$RNA, function(x) median(x$bcv_tagwise, na.rm = TRUE))
      corr_medians <- sapply(results$RNA_corr, function(x) median(x$bcv_tagwise, na.rm = TRUE))
      overall_red <- (median(rna_medians) - median(corr_medians)) / median(rna_medians) * 100
      cat(strrep("-", 60), "\n"); cat("COMPARISON (within-sample):\n")
      cat(sprintf("  Overall median BCV reduction: %.2f%%\n", overall_red))
      cat(sprintf("  samples showing reduction: %d/%d\n", sum(corr_medians < rna_medians), length(rna_medians)))
      if (overall_red > 0) cat("  -> RNA_corr shows REDUCED within-sample variation \n")
      else cat("  -> RNA_corr shows INCREASED within-sample variation \n")
    }
  }
  cat(strrep("=", 60), "\n\n")
}


#' Comprehensive BCV Analysis (Between- and Within-Sample)
#'
#' Performs a full Biological Coefficient of Variation (BCV) analysis for a list
#' of Seurat objects, including both between-sample and within-sample analyses.
#' Supports multiple correction methods for adjusted counts or offsets.
#'
#' @param objs Named list of \code{\link[Seurat]{Seurat}} objects.
#' @param output_dir Directory for saving outputs (plots, summaries, results).
#' @param min_counts_between Minimum counts threshold for between-sample BCV.
#' @param n_groups Number of cell groups per sample for within-sample BCV.
#' @param cells_per_group Optional fixed number of cells per group.
#' @param min_cells_per_group Minimum number of cells per group (within-sample).
#' @param min_counts_within Minimum counts threshold for within-sample BCV.
#' @param robust Logical; use robust dispersion estimation.
#' @param save_plots Logical; save generated plots to \code{output_dir}.
#' @param run_diagnostics Logical; run diagnostics to guide correction choice.
#' @param correction_method Correction method: \code{"auto"}, \code{"offset_ratio"},
#'   \code{"offset_od"}, \code{"direct"}, or \code{"offset_vector"}.
#' @param correction_vector_name Optional; feature metadata column to use for
#'   vector-based correction.
#' @param correction_vector_as Interpretation of \code{correction_vector_name}
#'   when using \code{"offset_vector"}: \code{"od"} (overdispersion values) or
#'   \code{"ratio"} (count ratios).
#' @param seed Random seed for reproducibility in within-sample grouping.
#'
#' @details
#' This function:
#' \enumerate{
#'   \item Runs optional diagnostics to determine the most appropriate correction method.
#'   \item Extracts pseudobulk counts for each assay.
#'   \item Calculates between-sample BCV for RNA and corrected counts (if available).
#'   \item Calculates within-sample BCV across cell groups per sample.
#'   \item Generates BCV plots for each analysis.
#'   \item Saves plots, RDS results, and summary tables if \code{save_plots = TRUE}.
#' }
#'
#' @return A list with components:
#' \describe{
#'   \item{results}{List with \code{between} and \code{within} BCV results.}
#'   \item{plots}{List of \code{ggplot2} objects generated.}
#'   \item{between_summary}{Summary table for between-sample BCV.}
#'   \item{within_summary}{Summary table for within-sample BCV.}
#'   \item{correction_method_used}{The correction method applied.}
#'   \item{correction_vector_name}{The name of the feature vector used for correction, if any.}
#'   \item{correction_vector_as}{Whether the correction vector was interpreted as "od" or "ratio".}
#' }
#'
#' @examples
#' \dontrun{
#' results <- analyze_bcv_comprehensive(
#'   objs = seurat_list,
#'   output_dir = "bcv_output",
#'   correction_method = "auto"
#' )
#' }
#'
#' @seealso
#' \code{\link{calculate_bcv_direct}}, \code{\link{analyze_within_sample_bcv}},
#' \code{\link{plot_bcv}}, \code{\link{plot_within_sample_bcv}},
#' \code{\link{plot_between_vs_within_bcv}}, \code{\link{create_summary_table}}
#'
#' @export


analyze_bcv_comprehensive <- function(
    objs,  # List of Seurat objects
    output_dir = "bcv_analysis_comprehensive",
    # Between-sample parameters
    min_counts_between = 10,
    # Within-sample parameters
    n_groups = 10,
    cells_per_group = NULL,
    min_cells_per_group = 50,
    min_counts_within = 10,
    # General parameters
    robust = TRUE,
    save_plots = TRUE,
    run_diagnostics = TRUE,
    correction_method = c("auto", "offset_ratio", "offset_od", "direct", "offset_vector"),
    correction_vector_name = NULL,
    correction_vector_as = c("od","ratio"),
    seed = 42
) {
  correction_method <- match.arg(correction_method)
  correction_vector_as <- match.arg(correction_vector_as)
  
  if (is.null(names(objs))) names(objs) <- paste0("sample_", seq_along(objs))
  if (save_plots && !dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  .log(strrep("=", 50))
  .log("Starting Comprehensive BCV Analysis")
  .log("Analyzing both between-sample and within-sample variation")
  .log(strrep("=", 50))
  
  results <- list(between = list(), within = list())
  plots <- list()
  
  # Diagnostics (optional)
  if (run_diagnostics) {
    .log("\n--- Running diagnostics ---")
    diag_results <- diagnose_seurat_correction(objs, sample_idx = 1)
    if (correction_method == "auto") {
      # Heuristic: if corrected/raw median ratio suggests RNA_corr is truly corrected,
      # prefer ratio offsets; otherwise direct.
      if (!is.null(diag_results)) {
        median_ratio <- median(diag_results$gene_ratios, na.rm = TRUE)
        if (median_ratio > 1.05) {
          .log("Auto-detected: corrected/raw > 1 -> using offset_ratio")
          correction_method <- "offset_ratio"
        } else {
          .log("Auto-detected: minimal/lower ratio -> using direct")
          correction_method <- "direct"
        }
      } else {
        correction_method <- "direct"
      }
    }
  }
  
  # Extract pseudobulks
  .log("\n", strrep("=", 40)); .log("BETWEEN-SAMPLE BCV ANALYSIS"); .log(strrep("=", 40))
  .log("\n--- Extracting pseudobulk data ---")
  pb_rna <- extract_and_pseudobulk(objs, assay_name = "RNA")
  has_rna_corr <- all(sapply(objs, function(obj) "RNA_corr" %in% names(obj@assays)))
  if (has_rna_corr) {
    pb_corr <- extract_and_pseudobulk(objs, assay_name = "RNA_corr")
    common_genes <- intersect(rownames(pb_rna), rownames(pb_corr))
    pb_rna <- pb_rna[common_genes, ]
    pb_corr <- pb_corr[common_genes, colnames(pb_rna)]
  } else {
    pb_corr <- NULL
  }
  
  # Optional feature vector (from feature_meta)
  feature_vec <- NULL
  if (!is.null(correction_vector_name)) {
    feature_vec <- extract_feature_vector(objs, correction_vector_name)
  }
  
  # RNA baseline
  .log("\n--- Calculating between-sample BCV for RNA (raw) ---")
  bcv_rna <- calculate_bcv_direct(pb_rna, min_counts = min_counts_between, robust = robust)
  results$between$RNA <- bcv_rna
  plots$between_RNA <- plot_bcv(bcv_rna, title = "Between-sample: RNA Assay", color = "darkblue")
  
  # Corrected path
  if (has_rna_corr || correction_method %in% c("offset_od","offset_vector")) {
    .log("\n--- Calculating between-sample BCV for corrected path ---")
    if (correction_method == "direct" && has_rna_corr) {
      .log("Using direct BCV from corrected counts (RNA_corr)")
      bcv_corr <- calculate_bcv_direct(pb_corr, min_counts = min_counts_between, robust = robust)
    } else if (correction_method == "offset_ratio" && has_rna_corr) {
      .log("Using offset correction based on count ratios (pb_corr/pb_rna)")
      bcv_corr <- calculate_bcv(pb_rna, pb_corr = pb_corr, offset_method = "ratio",
                                min_counts = min_counts_between, robust = robust)
    } else if (correction_method == "offset_od") {
      .log("Using offset correction based on overdispersion vector")
      od_values <- extract_overdispersion(objs)
      if (is.null(od_values)) stop("offset_od requested but no 'scaling_factor' found")
      bcv_corr <- calculate_bcv(pb_rna, overdispersion = od_values,
                                offset_method = "overdispersion",
                                min_counts = min_counts_between, robust = robust)
    } else if (correction_method == "offset_vector") {
      if (is.null(feature_vec)) stop("offset_vector requested but correction_vector_name not found")
      .log("Using offset from feature vector '%s' interpreted as '%s'",
           correction_vector_name, correction_vector_as)
      bcv_corr <- calculate_bcv(pb_rna, feature_vector = feature_vec,
                                vector_as = correction_vector_as,
                                offset_method = "vector",
                                min_counts = min_counts_between, robust = robust)
    } else {
      warning("Requested correction method not applicable; falling back to direct if RNA_corr present")
      if (has_rna_corr) {
        bcv_corr <- calculate_bcv_direct(pb_corr, min_counts = min_counts_between, robust = robust)
      } else {
        bcv_corr <- NULL
      }
    }
    
    if (!is.null(bcv_corr)) {
      results$between$RNA_corr <- bcv_corr
      plots$between_RNA_corr <- plot_bcv(bcv_corr, title = "Between-sample: Corrected", color = "darkgreen")
      plots$between_comparison <- plot_bcv_comparison(bcv_rna, bcv_corr)
    }
  }
  
  # ========== WITHIN-SAMPLE ==========
  .log("\n", strrep("=", 40)); .log("WITHIN-SAMPLE BCV ANALYSIS"); .log(strrep("=", 40))
  within_results <- analyze_within_sample_bcv(
    seurat_list = objs,
    assay_names = if (has_rna_corr) c("RNA","RNA_corr") else "RNA",
    n_groups = n_groups,
    cells_per_group = cells_per_group,
    min_cells_per_group = min_cells_per_group,
    min_counts = min_counts_within,
    robust = robust,
    seed = seed,
    offset_vector_name = correction_vector_name,  # NULL if not provided
    offset_vector_as = correction_vector_as
  )
  results$within <- within_results
  plots$within_summary <- plot_within_sample_summary(within_results)
  for (i in seq_along(names(objs))) {
    sample_name <- names(objs)[i]
    p <- plot_within_sample_bcv(within_results, sample_name = sample_name)
    if (!is.null(p)) plots[[paste0("within_", sample_name)]] <- p
  }
  
  # Combined
  if (!is.null(results$between$RNA_corr)) {
    plots$combined_comparison <- plot_between_vs_within_bcv(results$between, results$within)
  }
  
  # Save / summaries
  if (save_plots) {
    .log("\n--- Saving results ---")
    for (nm in grep("^between", names(plots), value = TRUE)) {
      ggsave(file.path(output_dir, paste0(nm, ".pdf")), plot = plots[[nm]], width = 8, height = 6); .log("Saved: %s.pdf", nm)
    }
    for (nm in grep("^within", names(plots), value = TRUE)) {
      ggsave(file.path(output_dir, paste0(nm, ".pdf")), plot = plots[[nm]],
             width = ifelse(grepl("summary", nm), 10, 12),
             height = ifelse(grepl("summary", nm), 8, 6))
      .log("Saved: %s.pdf", nm)
    }
    if ("combined_comparison" %in% names(plots)) {
      ggsave(file.path(output_dir, "combined_comparison.pdf"),
             plot = plots$combined_comparison, width = 10, height = 10)
      .log("Saved: combined_comparison.pdf")
    }
    saveRDS(results, file.path(output_dir, "bcv_results_comprehensive.rds"))
    between_summary <- create_summary_table(results$between)
    write.csv(between_summary, file.path(output_dir, "bcv_summary_between_sample.csv"), row.names = FALSE)
    within_summary <- summarize_within_sample_bcv(results$within)
    write.csv(within_summary, file.path(output_dir, "bcv_summary_within_sample.csv"), row.names = FALSE)
    .log("Saved summary tables")
  } else {
    between_summary <- create_summary_table(results$between)
    within_summary <- summarize_within_sample_bcv(results$within)
  }
  
  print_summary(results$between, type = "between")
  print_summary(results$within,  type = "within")
  
  cat("\n", strrep("=", 60), "\n"); cat("COMPREHENSIVE BCV ANALYSIS COMPLETE\n"); cat(strrep("=", 60), "\n")
  if (!is.null(results$between$RNA_corr)) {
    between_rna  <- median(results$between$RNA$bcv_tagwise, na.rm = TRUE)
    between_corr <- median(results$between$RNA_corr$bcv_tagwise, na.rm = TRUE)
    between_red  <- (between_rna - between_corr)/between_rna*100
    within_rna_medians  <- sapply(results$within$RNA,      function(x) median(x$bcv_tagwise, na.rm = TRUE))
    within_corr_medians <- if ("RNA_corr" %in% names(results$within))
      sapply(results$within$RNA_corr, function(x) median(x$bcv_tagwise, na.rm = TRUE)) else NA
    within_red <- if (all(!is.na(within_corr_medians)))
      (mean(within_rna_medians) - mean(within_corr_medians)) / mean(within_rna_medians) * 100 else NA
    cat("\nOVERALL IMPACT OF CORRECTION:\n")
    cat(sprintf("  Between-sample BCV reduction: %.1f%%\n", between_red))
    if (!is.na(within_red)) cat(sprintf("  Within-sample BCV reduction: %.1f%%\n", within_red))
    cat("\n")
  }
  cat(strrep("=", 60), "\n\n")
  
  list(
    results = results,
    plots = plots,
    between_summary = between_summary,
    within_summary = within_summary,
    correction_method_used = correction_method,
    correction_vector_name = correction_vector_name,
    correction_vector_as = correction_vector_as
  )
}