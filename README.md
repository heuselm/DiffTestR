# DiffTestR

An R package for testing differential protein abundance from AP-DIA-MS (Affinity Purification - Data-Independent Acquisition Mass Spectrometry) datasets.

## Overview

DiffTestR provides convenience functions for complete differential expression analysis workflows, from data import through statistical testing to visualization. The package handles precursor-level quantification data, performs robust normalization and imputation, and conducts statistical tests with protein-level rollup.

## Features

### Core Analysis
- **Flexible data import** from DIA-NN or custom table formats
- **Customizable normalization** (quantile, VSN, median, or custom methods)
- **Frequency-based filtering** to remove low-quality measurements
- **Smart missing value imputation** using low-abundance distributions
- **Statistical testing** at precursor level with protein-level summarization
- **Comprehensive visualizations** including volcano plots, heatmaps, and interactive HTML reports
- **Quality control plots** for normalization effects, sample correlations, and imputation diagnostics

### New in v0.9 (evolve branch)
- **Progress reporting** - Real-time feedback during long analyses with verbose mode
- **Summary statistics** - Automatic calculation and display of key metrics
- **S3 print methods** - Informative summary when printing results
- **Data quality assessment** - Comprehensive QC metrics with `assessDataQuality()`
  - Coefficient of variation analysis per condition
  - Leave-one-out CV analysis for outlier replicate detection
  - Missing value summaries
  - Dynamic range and replicate correlation metrics
- **Enhanced return objects** - Structured results with metadata and statistics
- **Multi-comparison testing** - NEW `testMultipleComparisons()` function
  - Auto-generates all pairwise comparisons from study design
  - OR accepts explicit comparison specification
  - Automated directory organization and result aggregation
  - Robust error handling across multiple comparisons
  - Combined visualization and summary statistics

## Installation

```r
# Install from GitHub (requires devtools)
devtools::install_github("heuselm/DiffTestR")

# Or install locally
devtools::install_local("path/to/DiffTestR")
```

## Dependencies

The package requires:
- R >= 3.50
- `data.table`, `ggplot2`, `corrplot`, `pheatmap`, `ggrepel`, `reshape2`
- `limma` (for normalization and statistics)
- `crosstalk`, `htmlwidgets`, `plotly` (for interactive visualizations)

## Quick Start

### 1. Import Data from DIA-NN

```r
library(DiffTestR)

# Import DIA-NN results (first run creates a study design template)
qdata <- importFromDIANN(
  path_to_result = "path/to/DIANN_results.tsv",
  quant_value_col = "Precursor.Quantity",
  study_design = "path/to/Study_design_filled.tsv"
)
```

**Study Design Format:**
Create a tab-separated file with columns:
- `filename`: Must match raw file names in your data
- `condition`: Experimental condition (e.g., "treated", "control")
- `replicate`: Replicate number (integer)

Example:
```
filename    condition    replicate
sample01    control      1
sample02    control      2
sample03    treated      1
sample04    treated      2
```

### 2. Test Differential Abundance

```r
# Run complete differential expression analysis
results <- testDifferentialAbundance(
  input_dt = qdata$data_wide,
  study_design = "path/to/Study_design_filled.tsv",
  condition_1 = "treated",
  condition_2 = "control",
  min_n_obs = 4,                          # Minimum observations per precursor
  normalize_data = TRUE,                   # Apply normalization
  normalization_function = limma::normalizeQuantiles,
  imp_percentile = 0.001,                 # Imputation centered at 0.1% quantile
  imp_sd = 0.2,                           # Imputation standard deviation
  plot_pdf = TRUE,                        # Generate diagnostic PDFs
  write_tsv_tables = TRUE,                # Write result tables
  target_protein = "Q9Y6K9"               # Optional: highlight specific protein
)
```

### 3. Run Multiple Comparisons (NEW in v0.9)

```r
# Auto-generate and test all pairwise comparisons
multi_result <- testMultipleComparisons(
  input_dt = qdata$data_wide,
  study_design = "Study_design.tsv",
  output_dir = "multi_comparison_results"  # Organized subdirectories
)

# View summary
print(multi_result)

# Access individual comparisons
multi_result$individual_results$ConditionA_vs_ConditionB

# Access combined results table
multi_result$combined_results_dt

# Or specify explicit comparisons
comparisons <- data.table(
  condition_1 = c("KO", "Treatment"),
  condition_2 = c("WT", "WT")
)

multi_result <- testMultipleComparisons(
  input_dt = qdata$data_wide,
  study_design = "Study_design.tsv",
  comparisons = comparisons
)
```

### 4. Visualize Multiple Comparisons

```r
# Create overview plots for multiple differential tests
overview <- plotDifferentialAbundanceOverview(
  diffExpr_result_dt = multi_result$combined_results_dt,
  significance_threshold_p = 0.05,
  significance_threshold_fc = 2,          # 2-fold change threshold
  browsable_html = TRUE,                   # Create interactive HTML plots
  heatmap = TRUE                          # Generate fold-change heatmap
)
```

## Detailed Workflow

### Step 1: Data Import

#### Option A: From DIA-NN Output

```r
qdata <- importFromDIANN(
  path_to_result = "DIANN_results.tsv",
  quant_value_col = "Precursor.Quantity",  # or "Precursor.Normalised"
  write_study_design_template = TRUE,
  study_design = "Study_design_filled.tsv"
)
```

The function returns a list containing:
- `data_long`: Long-format data table
- `data_wide`: Wide-format data table
- `matrix_raw`: Raw intensity matrix
- `ann_col`: Column annotations
- `ann_row`: Row annotations (protein/precursor info)

#### Option B: From Pre-formatted Tables

```r
qdata <- importFromTables(
  path_to_tsv = "precursor_matrix_wide.txt",
  precursor_annotation = "precursor_annotation.txt",
  study_design = "Study_design.txt"
)
```

### Step 2: Differential Expression Testing

The `testDifferentialAbundance` function performs:

1. **Data loading and filtering** to samples in study design
2. **Log2 transformation** of intensities
3. **Normalization** (default: quantile normalization)
4. **Quality control plots**: distributions, correlations, heatmaps
5. **Filtering** by minimum observations (`min_n_obs`)
6. **Missing value imputation** from low-abundance distribution
7. **T-tests** at precursor level
8. **Protein-level rollup** with Benjamini-Hochberg correction
9. **Volcano plots** at precursor and protein levels

```r
results <- testDifferentialAbundance(
  input_dt = qdata$data_wide,
  protein_group_annotation = NULL,        # Auto-extracted if NULL
  study_design = "Study_design.txt",
  normalize_data = TRUE,
  normalization_function = limma::normalizeQuantiles,  # Try also: normalizeVSN, normalizeMedianValues
  condition_1 = "treated",
  condition_2 = "control",
  min_n_obs = 4,
  imp_percentile = 0.001,
  imp_sd = 0.2,
  plot_pdf = TRUE,
  write_tsv_tables = TRUE,
  target_protein = "O08760"
)
```

**Output Files Generated:**
- `DiffTest_result.tsv`: Complete results table
- `Normalization_impact_abundance_density.pdf`: Normalization diagnostics
- `Heatmap_*.pdf`: Various heatmaps showing data quality
- `Correlation_sample_sample_*.pdf`: Sample correlation plots
- `Filtering_min_obs.pdf`: Distribution of observations
- `Imputation_*.pdf`: Imputation diagnostics
- `Volcano_plot_*.pdf`: Volcano plots at precursor and protein levels

**Results Object Contains:**
- `data_source`: Input file path
- `data_long`: Long-format processed data
- `data_matrix_log2`: Filtered, log2-transformed data
- `data_matrix_log2_imp`: Imputed data matrix
- `study_design`: Study design table
- `annotation_col`: Column annotations
- `diffExpr_result_dt`: Full results with statistics
- `candidates_condition1`: Significantly upregulated in condition 1
- `candidates_condition2`: Significantly upregulated in condition 2

### Step 3: Multi-Comparison Visualization

```r
# Combine results from multiple comparisons
all_results <- rbind(
  results1$diffExpr_result_dt,
  results2$diffExpr_result_dt,
  results3$diffExpr_result_dt
)

# Create overview plots
overview <- plotDifferentialAbundanceOverview(
  diffExpr_result_dt = all_results,
  significance_threshold_p = 0.05,
  significance_threshold_fc = 2,
  label_prefix = "MyExperiment",
  remove_from_protein_names = "\\(Bos;|_MOUSE|_HUMAN",
  browsable_html = TRUE,
  protein_highlight_tag = "IGSeq",
  heatmap = TRUE
)
```

**Output Files:**
- `MyExperiment.pdf`: Multi-panel volcano plots
- `MyExperiment_interactive_v1.html`: Basic interactive plot
- `MyExperiment_interactive_v2.html`: Interactive with cross-panel highlighting
- `MyExperiment_heatmap_diffProtFcProfiles.pdf`: Fold-change heatmap

## Advanced Options

### Custom Normalization

```r
# Use VSN normalization instead of quantile
results <- testDifferentialAbundance(
  input_dt = qdata$data_wide,
  study_design = "Study_design.txt",
  normalize_data = TRUE,
  normalization_function = limma::normalizeVSN
)

# Skip normalization (e.g., for pulldown vs input comparisons)
results <- testDifferentialAbundance(
  input_dt = qdata$data_wide,
  study_design = "Study_design.txt",
  normalize_data = FALSE
)
```

### Filtering and Imputation Parameters

```r
results <- testDifferentialAbundance(
  input_dt = qdata$data_wide,
  study_design = "Study_design.txt",
  min_n_obs = 6,              # More stringent filtering
  imp_percentile = 0.01,      # Center imputation at 1% quantile
  imp_sd = 0.3                # Wider imputation distribution
)
```

### Accessing Results

```r
# View top differentially abundant proteins
head(results$candidates_condition1)
head(results$candidates_condition2)

# Access full results table
diff_table <- results$diffExpr_result_dt

# Get specific proteins
target_results <- diff_table[Protein.Names %like% "MYC"]

# Export filtered candidates
write.table(results$candidates_condition1,
            "upregulated_proteins.tsv",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)
```

## Tips and Best Practices

1. **Replicates**: Include at least 3 biological replicates per condition for robust statistics
2. **Study Design**: Ensure filenames in study design exactly match those in your data
3. **Normalization**: For pulldown experiments comparing bait vs control, consider setting `normalize_data = FALSE`
4. **Filtering**: Adjust `min_n_obs` based on your total sample number and expected data completeness
5. **Imputation**: Lower `imp_percentile` values impute more conservatively (lower abundance)
6. **Quality Control**: Always review the generated PDFs to assess data quality and processing steps
7. **Reproducibility**: The analysis uses `set.seed(123)` for reproducible imputation

## Example Complete Analysis

```r
library(DiffTestR)

# 1. Import data
qdata <- importFromDIANN(
  path_to_result = "data/DIANN_results.tsv",
  study_design = "data/Study_design.tsv"
)

# 2. Test differential abundance
results <- testDifferentialAbundance(
  input_dt = qdata$data_wide,
  study_design = "data/Study_design.tsv",
  condition_1 = "KO",
  condition_2 = "WT",
  min_n_obs = 4,
  plot_pdf = TRUE
)

# 3. Examine results
print(paste("Found", nrow(results$candidates_condition1), "proteins upregulated in KO"))
print(paste("Found", nrow(results$candidates_condition2), "proteins upregulated in WT"))

# 4. Create overview (if multiple comparisons)
overview <- plotDifferentialAbundanceOverview(
  diffExpr_result_dt = results$diffExpr_result_dt,
  label_prefix = "KO_vs_WT"
)
```

## Citation

If you use DiffTestR in your research, please cite:

```
Heusel, M. DiffTestR: Testing Differential Expression from AP-DIA-MS datasets.
R package version 0.0.0.9000. https://github.com/heuselm/DiffTestR
```

## Author

Moritz Heusel (moritz.heusel@med.lu.se)
ORCID: 0000-0002-8506-530X

## License

MIT License
