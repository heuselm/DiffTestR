# Generate reference results for unit tests
# This script runs testDifferentialAbundance on the example data
# and saves the results as reference for regression testing

library(data.table)
library(limma)
library(ggplot2)
library(reshape2)
library(pheatmap)
library(corrplot)
library(ggrepel)
source("../../R/testDifferentialAbundance.R")

# Load example data
load("../../data/DiffTestR_example_data_wide.rda")
load("../../data/DiffTestR_example_study_design.rda")

# Run testDifferentialAbundance
set.seed(123) # For reproducible imputation
reference_results <- testDifferentialAbundance(
  input_dt = DiffTestR_example_data_wide,
  study_design = DiffTestR_example_study_design,
  condition_1 = "A",
  condition_2 = "B",
  min_n_obs = 4,
  normalize_data = TRUE,
  normalization_function = limma::normalizeQuantiles,
  imp_percentile = 0.001,
  imp_sd = 0.2,
  plot_pdf = FALSE,
  write_tsv_tables = FALSE
)

# Extract key statistics for regression testing
reference_stats <- list(
  n_precursors = nrow(reference_results$diffExpr_result_dt),
  n_proteins = length(unique(reference_results$diffExpr_result_dt$Protein.Group)),
  n_significant_precursors = sum(reference_results$diffExpr_result_dt$p_value_BHadj <= 0.05, na.rm = TRUE),
  n_significant_proteins = length(unique(reference_results$diffExpr_result_dt[p_value_BHadj_protein <= 0.05, Protein.Group])),
  median_log2fc_precursor = median(reference_results$diffExpr_result_dt$log2_fold_change, na.rm = TRUE),
  median_log2fc_protein = median(reference_results$diffExpr_result_dt$log2_fold_change_protein, na.rm = TRUE),
  # Sample a few specific protein results for exact matching
  sample_proteins = reference_results$diffExpr_result_dt[, .SD[1], by = Protein.Group][order(Protein.Group)][1:10]
)

# Save reference results
saveRDS(reference_stats, "reference_results.rds")

cat("Reference results generated:\n")
cat("  Precursors tested:", reference_stats$n_precursors, "\n")
cat("  Proteins tested:", reference_stats$n_proteins, "\n")
cat("  Significant precursors:", reference_stats$n_significant_precursors, "\n")
cat("  Significant proteins:", reference_stats$n_significant_proteins, "\n")
cat("  Median log2FC (precursor):", round(reference_stats$median_log2fc_precursor, 4), "\n")
cat("  Median log2FC (protein):", round(reference_stats$median_log2fc_protein, 4), "\n")
