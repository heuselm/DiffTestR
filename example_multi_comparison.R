# Example: Using testMultipleComparisons
# This script demonstrates how to use the testMultipleComparisons() function

library(DiffTestR)

# Load example data
data(DiffTestR_example_data_wide)
data(DiffTestR_example_study_design)

# Example 1: Auto-generate all pairwise comparisons
# =====================================================
# When comparisons = NULL, the function automatically generates
# all pairwise comparisons from unique conditions in the study design

result_auto <- testMultipleComparisons(
  input_dt = DiffTestR_example_data_wide,
  study_design = DiffTestR_example_study_design,
  comparisons = NULL,  # Auto-generate pairwise comparisons
  output_dir = "multi_comparison_auto",
  plot_pdf = FALSE,
  write_tsv_tables = TRUE,
  verbose = TRUE
)

# Print results summary
print(result_auto)

# Get detailed summary
summary(result_auto)

# Access individual results
names(result_auto$individual_results)

# Access combined results table
head(result_auto$combined_results_dt)


# Example 2: Explicit comparisons
# =====================================================
# You can also specify exactly which comparisons to perform

my_comparisons <- list(
  list(condition_1 = "A", condition_2 = "B")
)

result_explicit <- testMultipleComparisons(
  input_dt = DiffTestR_example_data_wide,
  study_design = DiffTestR_example_study_design,
  comparisons = my_comparisons,
  output_dir = "multi_comparison_explicit",
  plot_pdf = TRUE,
  write_tsv_tables = TRUE,
  normalize_data = TRUE,
  min_n_obs = 4,
  verbose = TRUE
)

print(result_explicit)


# Example 3: Access specific comparison results
# =====================================================

# Get results for a specific comparison
specific_result <- result_explicit$individual_results$A_vs_B

# This is a standard diffExpr object
print(specific_result)

# Access the result table
head(specific_result$diffExpr_result_dt)

# Get significant proteins
sig_proteins <- specific_result$diffExpr_result_dt[p_value_BHadj_protein <= 0.05]
nrow(sig_proteins)


# Example 4: Working with combined results
# =====================================================

# The combined_results_dt contains all results merged together
# with a 'comparison' column to identify which comparison each result belongs to

library(data.table)

# Count significant results per comparison
result_auto$combined_results_dt[
  p_value_BHadj_protein <= 0.05,
  .N,
  by = comparison
]

# Find proteins significant in multiple comparisons
sig_dt <- result_auto$combined_results_dt[p_value_BHadj_protein <= 0.05]
sig_proteins_by_comparison <- sig_dt[, .(comparisons = .N), by = Protein.Group]
proteins_sig_in_multiple <- sig_proteins_by_comparison[comparisons > 1]

# Get proteins with consistent direction across comparisons
protein_directions <- result_auto$combined_results_dt[
  p_value_BHadj_protein <= 0.05,
  .(
    mean_log2fc = mean(log2_fold_change_protein),
    n_comparisons = .N,
    all_same_direction = all(sign(log2_fold_change_protein) == sign(log2_fold_change_protein[1]))
  ),
  by = Protein.Group
]

# Show proteins significant and consistent across all comparisons
consistent_proteins <- protein_directions[all_same_direction == TRUE & n_comparisons > 1]
print(consistent_proteins)
