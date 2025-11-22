# DiffTestR Package Review & Improvement Suggestions

## Package Overview

**DiffTestR** is an R package for differential protein abundance testing from DIA-MS proteomics data. It implements a comprehensive workflow from data import through statistical testing with precursor-to-protein rollup.

### Current Strengths
1. ✅ **Complete workflow**: Import → Normalization → Filtering → Imputation → Testing → Rollup
2. ✅ **Flexible normalization**: Supports multiple normalization functions (quantile, VSN, median)
3. ✅ **Rich visualization**: Extensive QC plots, correlation matrices, heatmaps, volcano plots
4. ✅ **Interactive outputs**: HTML volcano plots with crosstalk for cross-panel highlighting
5. ✅ **Well-tested**: Comprehensive testthat suite with regression testing
6. ✅ **Example data**: Embedded 300-protein test dataset from HYE benchmark

---

## Critical Gaps in Information Flow

### 1. **Missing Progress Indicators** ⚠️ HIGH PRIORITY
**Problem**: Long-running analyses provide no progress feedback
- User has no idea if analysis is progressing or stuck
- No ETA for completion
- Difficult to debug when analyses fail mid-process

**Impact**: Poor UX, especially for large datasets (>5000 proteins)

**Suggested Solutions**:

#### Option A: Progress Bar with `progress` Package
```r
# Add to DESCRIPTION Imports: progress
library(progress)

# In testDifferentialAbundance():
pb <- progress_bar$new(
  format = "  Processing [:bar] :percent eta: :eta",
  total = nrow(data_to_process),
  clear = FALSE
)

# Update in loops:
for (i in seq_len(nrow(data))) {
  pb$tick()
  # ... processing
}
```

#### Option B: Verbose Logging with `cli` Package
```r
# Add to DESCRIPTION Imports: cli
library(cli)

cli_h1("DiffTestR Analysis Pipeline")
cli_alert_info("Loading data: {nrow(data)} precursors")
cli_progress_step("Normalizing data", spinner = TRUE)
# ... normalization code
cli_progress_done()
cli_alert_success("Normalization complete")
```

#### Option C: Simple Console Messages (Minimal, Already Available)
```r
# Add verbose parameter
testDifferentialAbundance <- function(..., verbose = TRUE) {
  if (verbose) {
    cat("\n[1/8] Loading data...\n")
    cat("  - Loaded", nrow(data), "precursors\n")
    cat("[2/8] Log2 transformation...\n")
    cat("[3/8] Quantile normalization...\n")
    # ... etc
  }
}
```

**Recommendation**: Start with **Option C** (simple), optionally upgrade to **Option B** (cli) for better formatting.

---

### 2. **No Summary Statistics Output** ⚠️ MEDIUM PRIORITY
**Problem**: Results tables lack summary metrics
- No quick overview of how many proteins/precursors tested
- No summary of significant hits
- No breakdown by fold-change magnitude

**Suggested Addition**:
```r
# Add to testDifferentialAbundance() return value:
summary_stats <- list(
  n_precursors_input = nrow(data),
  n_precursors_tested = nrow(result_dt),
  n_proteins_tested = length(unique(result_dt$Protein.Group)),
  n_significant_precursors = sum(result_dt$p_value_BHadj <= 0.05, na.rm = TRUE),
  n_significant_proteins = length(unique(result_dt[p_value_BHadj_protein <= 0.05, Protein.Group])),
  n_upregulated = sum(result_dt$log2_fold_change_protein > 1 &
                      result_dt$p_value_BHadj_protein <= 0.05, na.rm = TRUE),
  n_downregulated = sum(result_dt$log2_fold_change_protein < -1 &
                        result_dt$p_value_BHadj_protein <= 0.05, na.rm = TRUE),
  median_cv_condition1 = ...,
  median_cv_condition2 = ...
)

# Add print method:
print.diffExpr <- function(x, ...) {
  cat("\nDifferential Abundance Results\n")
  cat("===============================\n")
  cat("Precursors tested:", x$summary_stats$n_precursors_tested, "\n")
  cat("Proteins tested:", x$summary_stats$n_proteins_tested, "\n")
  cat("Significant proteins (p<0.05):", x$summary_stats$n_significant_proteins, "\n")
  cat("  Up-regulated:", x$summary_stats$n_upregulated, "\n")
  cat("  Down-regulated:", x$summary_stats$n_downregulated, "\n")
}
```

---

### 3. **Limited Data Quality Metrics** ⚠️ MEDIUM PRIORITY
**Problem**: No systematic quality assessment
- No CV calculation per condition
- No missing value summaries
- No dynamic range metrics
- No reproducibility metrics (between replicates)

**Suggested Addition**: New function `assessDataQuality()`
```r
#' Assess proteomics data quality
#'
#' @param data_matrix Matrix of log2 intensities
#' @param study_design Study design table
#' @return List with QC metrics
#' @export
assessDataQuality <- function(data_matrix, study_design) {

  # Calculate CVs per condition
  cv_by_condition <- calculate_cv_per_condition(data_matrix, study_design)

  # Missing value analysis
  missing_stats <- list(
    total_values = length(data_matrix),
    n_missing = sum(is.na(data_matrix)),
    pct_missing = 100 * sum(is.na(data_matrix)) / length(data_matrix),
    missing_per_sample = colSums(is.na(data_matrix)) / nrow(data_matrix) * 100
  )

  # Dynamic range
  dynamic_range <- list(
    min = min(data_matrix, na.rm = TRUE),
    max = max(data_matrix, na.rm = TRUE),
    range = max(data_matrix, na.rm = TRUE) - min(data_matrix, na.rm = TRUE)
  )

  # Replicate correlation
  replicate_cors <- calculate_replicate_correlations(data_matrix, study_design)

  return(list(
    cv_stats = cv_by_condition,
    missing_stats = missing_stats,
    dynamic_range = dynamic_range,
    replicate_correlations = replicate_cors
  ))
}
```

---

### 4. **No Intermediate Checkpoint Saving** ⚠️ LOW-MEDIUM PRIORITY
**Problem**: Long analyses can fail without recovery
- If analysis fails at testing step, all preprocessing lost
- No ability to resume from checkpoint
- Rerunning wastes computational time

**Suggested Addition**:
```r
testDifferentialAbundance <- function(..., checkpoint_dir = NULL) {

  if (!is.null(checkpoint_dir)) {
    dir.create(checkpoint_dir, showWarnings = FALSE)

    # Save after normalization
    saveRDS(normalized_data, file.path(checkpoint_dir, "01_normalized.rds"))

    # Save after filtering
    saveRDS(filtered_data, file.path(checkpoint_dir, "02_filtered.rds"))

    # Save after imputation
    saveRDS(imputed_data, file.path(checkpoint_dir, "03_imputed.rds"))

    # Final results
    saveRDS(results, file.path(checkpoint_dir, "04_results.rds"))
  }
}

# Recovery function
resumeAnalysis <- function(checkpoint_dir, from_step = "last") {
  # Find latest checkpoint and resume
}
```

---

### 5. **Limited Export Format Options** ⚠️ LOW PRIORITY
**Problem**: Only TSV export available
- No Excel export (common in biology)
- No JSON for integration with web tools
- No database-friendly format

**Suggested Addition**:
```r
#' Export results in multiple formats
#'
#' @param results DiffTestR results object
#' @param format One of "tsv", "xlsx", "json"
#' @param file Output file path
#' @export
exportResults <- function(results, format = "tsv", file) {
  switch(format,
    "tsv" = fwrite(results$diffExpr_result_dt, file, sep = "\t"),
    "xlsx" = writexl::write_xlsx(list(
      results = results$diffExpr_result_dt,
      summary = as.data.frame(results$summary_stats)
    ), file),
    "json" = jsonlite::write_json(results$diffExpr_result_dt, file, pretty = TRUE)
  )
}
```

---

### 6. **No Batch Effect Detection/Correction** ⚠️ MEDIUM PRIORITY (Advanced Feature)
**Problem**: No tools for batch effect assessment or correction
- Multi-batch experiments common in proteomics
- No PCA colored by batch
- No batch correction (e.g., ComBat)

**Suggested Addition**: New function `correctBatchEffects()`
```r
#' Detect and correct batch effects
#'
#' @param data_matrix Matrix of log2 intensities
#' @param study_design Study design with 'batch' column
#' @param method One of "combat", "limma", "none"
#' @return Batch-corrected matrix
#' @export
correctBatchEffects <- function(data_matrix, study_design, method = "combat") {

  if (!"batch" %in% names(study_design)) {
    stop("study_design must have 'batch' column")
  }

  switch(method,
    "combat" = {
      # Requires sva package
      require(sva)
      ComBat(dat = data_matrix,
             batch = study_design$batch,
             mod = model.matrix(~condition, data = study_design))
    },
    "limma" = {
      # Use limma::removeBatchEffect
      removeBatchEffect(data_matrix, batch = study_design$batch)
    },
    "none" = data_matrix
  )
}

# Visualization
plotBatchEffects <- function(data_matrix, study_design) {
  # PCA colored by batch vs condition
  # Before/after batch correction comparison
}
```

---

### 7. **Missing Power Analysis** ⚠️ LOW PRIORITY (Advanced Feature)
**Problem**: No guidance on sample size requirements
- Users don't know if they have enough replicates
- No power calculation for detecting fold changes

**Suggested Addition**:
```r
#' Estimate statistical power for differential abundance testing
#'
#' @param n_replicates Number of replicates per condition
#' @param fold_change Expected fold change to detect
#' @param cv Coefficient of variation (typical: 0.2-0.3)
#' @param alpha Significance level
#' @return Estimated power
#' @export
estimatePower <- function(n_replicates, fold_change, cv = 0.25, alpha = 0.05) {
  # Use power.t.test for simple case
  effect_size <- log2(fold_change) / cv

  power.t.test(
    n = n_replicates,
    delta = effect_size,
    sd = 1,
    sig.level = alpha,
    type = "two.sample"
  )
}
```

---

## Suggested New Features

### 8. **Multi-Comparison Support** ⚠️ MEDIUM-HIGH PRIORITY
**Current**: Only handles pairwise comparisons
**Improvement**: Support for multiple conditions in one run

```r
testDifferentialAbundance <- function(..., comparisons = NULL) {

  # comparisons = list(
  #   c("A", "B"),
  #   c("A", "C"),
  #   c("B", "C")
  # )

  if (!is.null(comparisons)) {
    results_list <- lapply(comparisons, function(comp) {
      testDifferentialAbundance(
        ...,
        condition_1 = comp[1],
        condition_2 = comp[2]
      )
    })

    # Combine results
    combined_results <- combineComparisons(results_list)
    return(combined_results)
  }
}
```

### 9. **Integration with Pathway Analysis** ⚠️ LOW PRIORITY (Future)
**Idea**: Direct integration with enrichment tools
```r
#' Perform pathway enrichment on significant proteins
#'
#' @param difftest_results Results from testDifferentialAbundance
#' @param database One of "GO", "KEGG", "Reactome"
#' @export
enrichPathways <- function(difftest_results, database = "GO") {
  # Extract significant proteins
  sig_proteins <- difftest_results$diffExpr_result_dt[
    p_value_BHadj_protein <= 0.05,
    unique(Protein.Group)
  ]

  # Interface with clusterProfiler or similar
  # Return enrichment results
}
```

---

## Implementation Priority Ranking

### **Must Have** (Immediate Impact)
1. **Progress indicators** (#1) - Dramatically improves UX
2. **Summary statistics** (#2) - Essential for result interpretation
3. **Multi-comparison support** (#8) - Common user need

### **Should Have** (High Value)
4. **Data quality metrics** (#3) - Important for method validation
5. **Batch effect detection** (#6) - Common in real-world data

### **Nice to Have** (Lower Priority)
6. **Checkpoint saving** (#4) - Helpful for large datasets
7. **Export format options** (#5) - Convenience feature
8. **Power analysis** (#7) - Advanced users
9. **Pathway integration** (#9) - Future enhancement

---

## Architectural Suggestions

### 1. **Modularize the Pipeline**
Current `testDifferentialAbundance()` is monolithic (~600 lines). Consider breaking into:
```r
testDifferentialAbundance <- function(...) {
  # Orchestrator function
  data <- loadAndValidateData(...)
  data_norm <- normalizeData(data, ...)
  data_filt <- filterByObservations(data_norm, ...)
  data_imp <- imputeMissingValues(data_filt, ...)
  results <- performStatisticalTests(data_imp, ...)
  results <- rollupToProteinLevel(results, ...)

  if (plot_pdf) generateQCPlots(data_norm, data_filt, data_imp)

  return(results)
}
```

**Benefits**:
- Easier to test individual steps
- Users can call sub-functions for custom workflows
- Better code maintainability

### 2. **Add S3 Class Methods**
```r
# Create diffExpr S3 class
class(results) <- c("diffExpr", "list")

# Methods
print.diffExpr()
summary.diffExpr()
plot.diffExpr()
```

### 3. **Consistent Return Structure**
Ensure all functions return standardized objects with:
- **data**: processed data
- **metadata**: analysis parameters
- **stats**: summary statistics
- **qc**: quality metrics

---

## Discussion Points for @rcoder

### Questions for Implementation:

1. **Progress Reporting**: Which approach? Simple cat() messages or invest in `cli` package?

2. **Backward Compatibility**: Should we maintain exact current behavior with opt-in for new features?

3. **Dependencies**: Willing to add new dependencies (cli, progress, sva for batch correction)?

4. **API Breaking Changes**: Some suggestions require function signature changes. Version bump to 1.0.0?

5. **Performance**: Any computationally intensive steps that need optimization (parallelization)?

6. **Documentation**: Should we create vignettes for:
   - Basic workflow
   - Handling batch effects
   - Multi-condition comparisons
   - Quality control best practices

7. **Validation**: How to validate new features? More test data needed?

---

## Immediate Action Items

If implementing improvements, I recommend this order:

### Phase 1 (Quick Wins - 1-2 days)
- [ ] Add verbose parameter with console messages (#1 Option C)
- [ ] Add summary statistics to return object (#2)
- [ ] Add print.diffExpr() method
- [ ] Update documentation

### Phase 2 (Core Improvements - 3-5 days)
- [ ] Implement assessDataQuality() function (#3)
- [ ] Add multi-comparison support (#8)
- [ ] Modularize testDifferentialAbundance() pipeline
- [ ] Write vignette for basic workflow

### Phase 3 (Advanced Features - 1 week+)
- [ ] Batch effect detection/correction (#6)
- [ ] Checkpoint/resume functionality (#4)
- [ ] Additional export formats (#5)
- [ ] Power analysis tools (#7)

---

## Conclusion

DiffTestR is a solid, functional package with excellent statistical methodology. The main gaps are in **user experience** and **real-world workflow support**:

- Users need **feedback** during long analyses
- Users need **quality metrics** to trust results
- Users need **flexibility** for complex experimental designs

Implementing the high-priority items (#1, #2, #3, #8) would significantly enhance the package's usability without major architectural changes.

**Ready to discuss priorities with @rcoder!**
