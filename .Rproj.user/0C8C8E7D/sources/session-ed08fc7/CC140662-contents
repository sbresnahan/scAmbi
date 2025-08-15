#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List compute_overdisp_block_improved(S4 boot_sparse,
                                     int start_col, int end_col,
                                     int n_boot, int n_genes,
                                     double eps = 1e-8,
                                     double pseudocount = 0.1) {
  IntegerVector i = boot_sparse.slot("i");  // row indices
  IntegerVector p = boot_sparse.slot("p");  // column pointers (len ncol+1)
  NumericVector x = boot_sparse.slot("x");  // nonzero values
  
  const int ncol_total = p.size() - 1; // number of columns
  const int nnz        = x.size();
  
  if (n_boot <= 1) Rcpp::stop("n_boot must be >= 2.");
  if (start_col < 0) start_col = 0;
  if (end_col > ncol_total) end_col = ncol_total; // EXCLUSIVE
  if (start_col >= end_col)
    Rcpp::stop("start_col must be < end_col (exclusive) within [0, ncol].");
  
  if (p[0] != 0) Rcpp::stop("Invalid column pointer: p[0] != 0.");
  if (p[ncol_total] != nnz)
    Rcpp::stop(std::string("Invalid column pointer: p[last]=")
                 + std::to_string(p[ncol_total]) + " != length(x)="
                 + std::to_string(nnz) + ".");
                 
                 const int span = end_col - start_col;
                 const int bs   = span / n_boot; // full cells only in this block
                 if (bs <= 0) Rcpp::stop("No complete bootstrap blocks in the requested range.");
                 
                 NumericVector OD_sum(n_genes, 0.0);
                 IntegerVector DF_sum(n_genes, 0);
                 NumericVector expr_cells(n_genes, 0.0);
                 
                 // Assume cell-major layout: [cell0: b0..bN-1 | cell1: b0..bN-1 | ...]
                 for (int cell = 0; cell < bs; ++cell) {
                   const int base_col = start_col + cell * n_boot;
                   
                   // Accumulate across the n_boot columns for this cell.
                   std::vector<double> gene_sums(n_genes, 0.0);
                   std::vector<double> gene_sums_sq(n_genes, 0.0);
                   std::vector<char>   gene_seen(n_genes, 0); // seen >= once across boots
                   
                   for (int b = 0; b < n_boot; ++b) {
                     const int col = base_col + b;
                     if (col + 1 > end_col) break;
                     
                     const int a   = p[col];
                     const int bnd = p[col + 1];
                     if (a < 0 || bnd < a || bnd > nnz) {
                       Rcpp::stop(std::string("Invalid p-range for col=") + std::to_string(col)
                                    + ": p[col]=" + std::to_string(a)
                                    + ", p[col+1]=" + std::to_string(bnd)
                                    + ", nnz=" + std::to_string(nnz));
                     }
                     
                     for (int idx = a; idx < bnd; ++idx) {
                       const int g = i[idx];
                       if ((unsigned)g >= (unsigned)n_genes) continue;
                       const double val = x[idx];
                       if (val != 0.0) gene_seen[g] = 1; // zeros are implicit in sparse format
                       gene_sums[g]    += val;
                       gene_sums_sq[g] += val * val;
                     }
                   }
                   
                   // Convert sums to per-gene OD for this cell (include zeros via k=n_boot)
                   for (int g = 0; g < n_genes; ++g) {
                     if (!gene_seen[g]) continue;     // entirely zero across boots for this cell
                     const int    k   = n_boot;       // include zeros in sample size
                     const double mu  = gene_sums[g] / (double) k;
                     // Unbiased sample variance across the n_boot replicates
                     const double var = (gene_sums_sq[g] - (gene_sums[g] * gene_sums[g]) / (double) k)
                       / (double)(k - 1);
                     
                     // Inflation factor: OD = 1 + (bootstrap variance)/(mean)
                     const double od = 1.0 + (var / (mu + eps));
                     
                     OD_sum[g]     += od;   // one OD observation from this cell
                     DF_sum[g]     += 1;    // count contributing cells
                     expr_cells[g] += 1.0;  // same as DF_sum but numeric
                   }
                 }
                 
                 return List::create(
                   Named("OD_sum")     = OD_sum,
                   Named("DF_sum")     = DF_sum,
                   Named("expr_cells") = expr_cells
                 );
}