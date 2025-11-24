# =============================================================================
# DiffTestR ddata S3 Class - Usage Examples
# =============================================================================
#
# This script demonstrates the usage of the ddata S3 class methods
# for working with differential expression results from testDifferentialAbundance()
#
# =============================================================================

library(DiffTestR)

# =============================================================================
# 1. Basic Workflow: Import -> Test -> Explore
# =============================================================================

# Step 1: Import data to create qdata object
qdata <- importFromTables(
  path_to_tsv = "data/precursor_matrix.tsv",
  study_design = "data/study_design.txt"
)

# Step 2: Perform differential abundance testing
# This returns a ddata object
ddata <- testDifferentialAbundance(
  data = qdata,
  condition_1 = "Control",
  condition_2 = "Treatment",
  normalize_data = TRUE,
  min_n_obs = 4,
  verbose = TRUE
)

# =============================================================================
# 2. Print Method - Concise Summary
# =============================================================================

# Simply printing the object shows a concise summary
print(ddata)

# Or just type the object name
ddata

# =============================================================================
# 3. Summary Method - Detailed Statistics
# =============================================================================

# Get detailed summary with top candidate proteins
summary(ddata)

# =============================================================================
# 4. Plot Method - Volcano Plots
# =============================================================================

# Default: Static volcano plot with ggplot2
plot(ddata)

# Interactive volcano plot with plotly
plot(ddata, interactive = TRUE)

# Customize thresholds
plot(ddata, p_cutoff = 0.01, fc_cutoff = 2)

# No labels
plot(ddata, label_top = 0)

# More labels
plot(ddata, label_top = 20)

# Adjust point appearance
plot(ddata, point_size = 3, alpha = 0.8)

# =============================================================================
# 5. Plot Method - MA Plots
# =============================================================================

# Static MA plot
plot(ddata, which = "ma")

# Interactive MA plot
plot(ddata, which = "ma", interactive = TRUE)

# Customize MA plot
plot(ddata, which = "ma", p_cutoff = 0.01, fc_cutoff = 1.5, label_top = 15)

# =============================================================================
# 6. Accessing ddata Components
# =============================================================================

# Access the full results table
head(ddata$diffExpr_result_dt)

# Get candidate proteins higher in condition 1
candidates_cond1 <- ddata$candidates_condition1
print(candidates_cond1)

# Get candidate proteins higher in condition 2
candidates_cond2 <- ddata$candidates_condition2
print(candidates_cond2)

# Access summary statistics
ddata$summary_stats

# Access comparison information
ddata$comparison
ddata$conditions$condition_1
ddata$conditions$condition_2

# Access processed matrices
dim(ddata$mat_quant_log2_qnorm_imp_minObs)
head(ddata$mat_quant_log2_qnorm_imp_minObs[, 1:3])

# Study design
ddata$study_design

# =============================================================================
# 7. Filtering and Subsetting Results
# =============================================================================

library(data.table)

# Get all significant proteins (p < 0.05)
sig_proteins <- unique(ddata$diffExpr_result_dt[
  p_value_BHadj_protein <= 0.05,
  .(Protein.Group, Protein.Names, log2_fold_change_protein, p_value_BHadj_protein)
])
print(sig_proteins)

# Get highly regulated proteins (|log2FC| > 2, p < 0.01)
highly_regulated <- unique(ddata$diffExpr_result_dt[
  p_value_BHadj_protein <= 0.01 & abs(log2_fold_change_protein) > 2,
  .(Protein.Group, Protein.Names, log2_fold_change_protein, p_value_BHadj_protein)
])
print(highly_regulated)

# Search for specific proteins
myc_proteins <- ddata$diffExpr_result_dt[
  grepl("MYC", Protein.Names, ignore.case = TRUE)
]
print(myc_proteins)

# =============================================================================
# 8. Validate ddata Structure
# =============================================================================

# Validate that the object has all required components
validate_ddata(ddata)

# This is useful for checking objects after loading from file
# or when creating custom ddata objects

# =============================================================================
# 9. Saving and Loading Results
# =============================================================================

# Save the entire ddata object
saveRDS(ddata, "results/ddata_results.rds")

# Load the ddata object
ddata_loaded <- readRDS("results/ddata_results.rds")

# Verify it's still valid
validate_ddata(ddata_loaded)

# Use the loaded object
print(ddata_loaded)
plot(ddata_loaded, interactive = TRUE)

# =============================================================================
# 10. Export Results to Files
# =============================================================================

# Export full results table to TSV
write.table(
  ddata$diffExpr_result_dt,
  "results/differential_expression_results.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# Export significant proteins only
sig_results <- unique(ddata$diffExpr_result_dt[
  p_value_BHadj_protein <= 0.05,
  .(Protein.Group, Protein.Names, log2_fold_change_protein,
    p_value_BHadj_protein, n_precursors)
])

write.table(
  sig_results,
  "results/significant_proteins.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# Export candidates for each condition
write.table(
  ddata$candidates_condition1,
  sprintf("results/candidates_higher_in_%s.tsv", ddata$conditions$condition_1),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  ddata$candidates_condition2,
  sprintf("results/candidates_higher_in_%s.tsv", ddata$conditions$condition_2),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# =============================================================================
# 11. Combining with Other Plots
# =============================================================================

# Create a custom volcano plot using the data
library(ggplot2)

# Extract protein-level data
plot_data <- unique(ddata$diffExpr_result_dt[, .(
  Protein.Group,
  Protein.Names,
  log2_fold_change_protein,
  p_value_BHadj_protein
)])

# Custom ggplot
custom_volcano <- ggplot(plot_data, aes(
  x = log2_fold_change_protein,
  y = -log10(p_value_BHadj_protein)
)) +
  geom_point(aes(color = p_value_BHadj_protein <= 0.05 & abs(log2_fold_change_protein) >= 1)) +
  scale_color_manual(values = c("grey", "red")) +
  theme_minimal() +
  labs(title = paste("Custom Volcano:", ddata$comparison))

print(custom_volcano)

# =============================================================================
# 12. Working with Multiple Comparisons
# =============================================================================

# If you have multiple ddata objects from different comparisons,
# you can combine results

# Example: Three comparisons
ddata1 <- testDifferentialAbundance(qdata, condition_1 = "A", condition_2 = "B")
ddata2 <- testDifferentialAbundance(qdata, condition_1 = "A", condition_2 = "C")
ddata3 <- testDifferentialAbundance(qdata, condition_1 = "B", condition_2 = "C")

# Combine significant proteins across comparisons
all_sig <- rbind(
  ddata1$diffExpr_result_dt[p_value_BHadj_protein <= 0.05, .(Protein.Group, comparison)],
  ddata2$diffExpr_result_dt[p_value_BHadj_protein <= 0.05, .(Protein.Group, comparison)],
  ddata3$diffExpr_result_dt[p_value_BHadj_protein <= 0.05, .(Protein.Group, comparison)]
)

# Find proteins significant in all comparisons
sig_in_all <- all_sig[, .N, by = Protein.Group][N == 3]
print(sig_in_all)

# =============================================================================
# End of Examples
# =============================================================================
