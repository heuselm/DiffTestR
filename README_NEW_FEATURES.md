# New Features in DiffTestR v0.9 (evolve branch)

## Progress Reporting

All analysis functions now support verbose progress reporting to provide real-time feedback during long-running analyses.

```r
# Run analysis with progress messages (default)
results <- testDifferentialAbundance(
  input_dt = my_data,
  study_design = my_design,
  condition_1 = "A",
  condition_2 = "B",
  verbose = TRUE  # default
)

# Output:
# ============================================================
#  DiffTestR: Differential Abundance Analysis
# ============================================================
#
# [1/8] Loading data
#   Loaded 15234 precursors from R object
# [2/8] Log2 transformation
# [3/8] Quantile normalization
# [4/8] Filtering (min_n_obs >= 4)
#   Retained 12456/15234 precursors after filtering
# [5/8] Imputing missing values
#   Imputed 3421 missing values
# [6/8] Precursor-level statistical testing (12456 tests)
# [7/8] Protein-level rollup and multiple testing correction
#   Found 342 significant proteins (p < 0.05)
# [8/8] Finalizing results
# [OK] Analysis complete!
```

Disable progress messages for non-interactive use:
```r
results <- testDifferentialAbundance(..., verbose = FALSE)
```

## Enhanced Result Objects

Results now include comprehensive summary statistics and use S3 classes for better printing:

```r
results <- testDifferentialAbundance(...)

# Print automatically shows summary
print(results)

# Output:
# ====================================
#   Differential Abundance Results
# ====================================
#
# Input Data:
#   Precursors (input): 15234
#   Precursors (tested): 12456
#   Proteins (tested): 2341
#
# Significant Results (p < 0.05):
#   Precursors: 1834 (14.7%)
#   Proteins: 342 (14.6%)
#
# Regulation (|log2FC| > 1, p < 0.05):
#   Up-regulated: 156 proteins
#   Down-regulated: 142 proteins
#
# Fold Changes:
#   Median log2FC (protein): 0.234
#
# Data Quality (CV %):
#   Condition 1: 18.3%
#   Condition 2: 21.1%
#
# Comparison:
#   A vs B
# ====================================
```

Access summary statistics programmatically:
```r
results$summary_stats$n_significant_proteins
# [1] 342

results$summary_stats$median_cv_condition1
# [1] 18.3
```

## Data Quality Assessment

New comprehensive quality control function with leave-one-out CV analysis for outlier detection:

```r
# Prepare data matrix
data_mat <- as.matrix(my_data[, 4:11])
rownames(data_mat) <- my_data$Precursor.Id
data_mat_log2 <- log2(data_mat + 1)

# Run comprehensive QC
qc <- assessDataQuality(
  data_matrix = data_mat_log2,
  study_design = my_study_design,
  cv_threshold = 0.3,  # Flag replicates if removing improves CV by >30%
  verbose = TRUE
)

# Output:
# ============================================================
#  Data Quality Assessment
# ============================================================
#
# --- Coefficient of Variation Analysis ---
# [INFO] Condition 'A': Median CV = 18.3%
# [INFO] Condition 'B': Median CV = 21.1%
#
# --- Leave-One-Out Replicate Analysis ---
# [WARNING] Outlier detected: sample_A_rep3 (CV improves by 42.1% when removed)
# [OK] Condition 'B': No outlier replicates detected
#
# --- Missing Value Analysis ---
# [INFO] Overall missing values: 12.4%
#
# --- Dynamic Range ---
# [INFO] Dynamic range: 8.23 to 25.67 (range: 17.44 log2 units)
#
# --- Replicate Correlation ---
# [INFO] Condition 'A': Median correlation = 0.945
# [INFO] Condition 'B': Median correlation = 0.963
```

### QC Results Structure

```r
# Access CV statistics
qc$cv_stats$A$median_cv
# [1] 18.3

# Check for outlier replicates
print(qc$outlier_replicates)
#    condition              filename replicate baseline_cv    loo_cv cv_improvement is_outlier
# 1:         A sample_A_rep3         3       18.3      10.6          0.421       TRUE

# Missing value statistics
qc$missing_stats$pct_missing
# [1] 12.4

# Replicate correlations
qc$replicate_correlations$A$median_cor
# [1] 0.945
```

### Leave-One-Out CV Analysis

The LOO CV analysis helps identify problematic replicates:

```r
# For each condition with ≥3 replicates:
# 1. Calculate baseline CV across all replicates
# 2. Recalculate CV after removing each replicate individually
# 3. Flag replicate if its removal improves median CV by >threshold

# Example interpretation:
# If sample_A_rep3 is flagged with 42% improvement:
# - Baseline CV (all 4 reps): 18.3%
# - CV without rep3 (3 reps): 10.6%
# - Improvement: (18.3 - 10.6) / 18.3 = 42.1%
# → Consider excluding rep3 from analysis
```

## Complete Example Workflow

```r
library(DiffTestR)

# 1. Load example data
data(DiffTestR_example_data_wide)
data(DiffTestR_example_study_design)

# 2. Quality assessment before analysis
data_mat <- as.matrix(DiffTestR_example_data_wide[, 4:11])
rownames(data_mat) <- DiffTestR_example_data_wide$Precursor.Id
data_mat_log2 <- log2(data_mat + 1)

qc_pre <- assessDataQuality(data_mat_log2, DiffTestR_example_study_design)
print(qc_pre)

# Check for outliers
if (nrow(qc_pre$outlier_replicates) > 0) {
  cat("\nWarning: Outlier replicates detected:\n")
  print(qc_pre$outlier_replicates[, .(condition, filename, cv_improvement)])
}

# 3. Run differential abundance analysis with progress reporting
results <- testDifferentialAbundance(
  input_dt = DiffTestR_example_data_wide,
  study_design = DiffTestR_example_study_design,
  condition_1 = "A",
  condition_2 = "B",
  min_n_obs = 4,
  normalize_data = TRUE,
  plot_pdf = TRUE,
  write_tsv_tables = TRUE,
  verbose = TRUE  # Show progress
)

# 4. View results summary
print(results)

# 5. Access specific results
sig_proteins <- results$diffExpr_result_dt[p_value_BHadj_protein <= 0.05]
cat("\nTop 10 significant proteins:\n")
print(sig_proteins[order(p_value_BHadj_protein)][1:10,
                   .(Protein.Group, log2_fold_change_protein,
                     p_value_BHadj_protein, n_precursors)])

# 6. Extract candidates
up_regulated <- results$candidates_condition1
down_regulated <- results$candidates_condition2

cat(sprintf("\nUp-regulated in %s: %d proteins\n",
            results$summary_stats$condition_1,
            nrow(up_regulated)))
cat(sprintf("Down-regulated in %s: %d proteins\n",
            results$summary_stats$condition_1,
            nrow(down_regulated)))
```

## API Changes

### Backward Compatibility

All new features are **backward compatible**. Existing code will continue to work without modification.

### New Parameters

**`testDifferentialAbundance()`**:
- `verbose = TRUE` - Control progress reporting (new, default TRUE)

**`assessDataQuality()`** (new function):
- `data_matrix` - Matrix of log2 intensities
- `study_design` - Study design data.table
- `cv_threshold = 0.3` - Threshold for flagging outlier replicates
- `verbose = TRUE` - Control progress messages

### Enhanced Return Objects

**`testDifferentialAbundance()` now returns**:
- All previous components (unchanged)
- `summary_stats` - List of summary statistics (new)
- S3 class `"diffExpr"` for enhanced printing (new)

**`assessDataQuality()` returns**:
- `cv_stats` - CV statistics per condition
- `cv_loo_analysis` - Leave-one-out CV analysis results
- `outlier_replicates` - Data.table of flagged outliers
- `missing_stats` - Missing value analysis
- `dynamic_range` - Min, max, range of intensities
- `replicate_correlations` - Correlation matrices per condition
- S3 class `"DiffTestR_QC"` for enhanced printing

## Performance Notes

- Progress reporting adds negligible overhead (<0.1% runtime)
- Summary statistics calculated only once at end of analysis
- QC assessment scales linearly with number of features and samples
- Leave-one-out CV analysis requires ≥3 replicates per condition

## Migration Guide

No changes required for existing code! New features are opt-in or automatic enhancements.

To adopt new features:
```r
# Before (still works)
results <- testDifferentialAbundance(input_dt, study_design, "A", "B")

# After (enhanced experience)
results <- testDifferentialAbundance(input_dt, study_design, "A", "B", verbose = TRUE)
print(results)  # Now shows summary
qc <- assessDataQuality(log2_matrix, study_design)  # New QC function
```
