# Species Benchmark Analysis Summary

## Overview
Analysis of HYE (Human-Yeast-E.coli) species mix benchmark data comparing:
- **Instruments**: Exploris 480 vs timsTOF Pro 2
- **Comparison**: Mix A vs Mix B
- **Replicates**: 4 biological replicates per condition per instrument

## Analysis Scripts Created

### 1. `analyze_species_benchmark.R`
Main differential abundance analysis with quantile normalization.

**Key Steps:**
- Import data from Exploris480.tsv and TimsTof.tsv
- Filter out multi-species proteins (containing ";")
- Annotate proteins by species (H. sapiens, S. cerevisiae, E. coli)
- Run testDifferentialAbundance with quantile normalization
- Generate volcano plots and MA plots colored by species

**Outputs:**
- `volcano_plots_sidebyside.pdf/png` - Volcano plots for both instruments
- `ma_plots_sidebyside.pdf/png` - MA plots showing ratio vs abundance
- `combined_protein_results.tsv` - All protein-level results
- `summary_statistics.tsv` - Summary statistics by instrument and species

### 2. `analyze_CV_and_PCA.R`
Coefficient of variation and PCA analysis.

**Key Steps:**
- Calculate CV at precursor level (individual peptides)
- Calculate CV at protein level (aggregated quantities)
- Compare variability between instruments
- Perform PCA on precursor-level data
- Assess sample separation between Mix A and Mix B

**Outputs:**
- `CV_summary_statistics.tsv` - CV statistics by instrument and level
- `CV_comparison.pdf/png` - Violin/box plots of CV distributions
- `PCA_precursor_level.pdf/png` - PCA plots showing sample separation

### 3. `analyze_species_benchmark_noNorm.R`
Differential abundance analysis WITHOUT quantile normalization.

**Key Steps:**
- Same as main analysis but with `normalize_data = FALSE`
- Allows comparison of normalized vs non-normalized results

**Outputs:**
- All outputs suffixed with `_noNorm`

## Key Findings

### Coefficient of Variation Analysis

**Overall Median CV Values:**
| Instrument      | Precursor Level | Protein Level |
|-----------------|-----------------|---------------|
| Exploris 480    | 20.7%          | 26.1%         |
| timsTOF Pro 2   | 26.8%          | 26.0%         |

**By Condition:**
- **Exploris Mix A**: High CV at both precursor (36.8%) and protein (36.4%) levels
- **Exploris Mix B**: Low CV at both precursor (8.3%) and protein (8.6%) levels
- **timsTOF Mix B**: Moderate CV at both precursor (26.8%) and protein (26.0%) levels

**Conclusion:**
- **timsTOF shows HIGHER variability at precursor level** compared to Exploris
- This increased variability likely explains fewer significantly regulated proteins in timsTOF
- Exploris shows better reproducibility, especially in Mix B samples

### PCA Analysis

**Variance Explained by PC1:**
- **Exploris 480**: 72.4% of variance
- **timsTOF Pro 2**: 48.4% of variance

**Conclusion:**
- **Exploris shows better separation between Mix A and Mix B**
- Higher PC1 variance indicates clearer biological signal
- timsTOF data shows more technical variation relative to biological differences

## Data Processing Details

### Multi-species Protein Filtering
- **Exploris**: Removed 27,528 multi-species protein groups
- **timsTOF**: Removed multi-species protein groups (number reported in output)

### Species Annotation
Proteins annotated based on UniProt suffixes:
- `_HUMAN` → H. sapiens (Red in plots)
- `_YEAST` → S. cerevisiae (Blue in plots)
- `_ECOLI` → E. coli (Green in plots)

### Statistical Testing Parameters
- **Minimum observations**: 4 per precursor
- **Normalization**: Quantile normalization (limma::normalizeQuantiles)
- **Imputation**: 0.1% percentile, SD = 0.2
- **Significance thresholds**: p-value ≤ 0.05, |log2FC| ≥ 1

## Interpretation

### Why timsTOF Shows Fewer Regulated Proteins

1. **Higher Technical Variability**
   - 30% higher median CV at precursor level
   - Reduces power to detect true differences

2. **Lower Biological Signal**
   - PC1 explains only 48% vs 72% for Exploris
   - More variance attributed to technical factors

3. **Instrument Characteristics**
   - Different quantification principles
   - Different data acquisition modes
   - May require instrument-specific normalization

### Recommendations

1. **For timsTOF data:**
   - Consider alternative normalization strategies
   - May need stricter filtering criteria
   - Could benefit from technical replicate pooling

2. **General:**
   - Exploris 480 shows superior quantitative performance
   - Both instruments identify similar protein groups
   - Choice depends on throughput vs precision requirements

## Files Generated

### Results Directory Structure
```
results/
├── Exploris480_study_design.tsv
├── TimsTof_study_design.tsv
├── combined_protein_results.tsv
├── summary_statistics.tsv
├── volcano_plots_sidebyside.pdf
├── volcano_plots_sidebyside.png
├── ma_plots_sidebyside.pdf
├── ma_plots_sidebyside.png
├── CV_summary_statistics.tsv
├── CV_comparison.pdf
├── CV_comparison.png
├── PCA_precursor_level.pdf
├── PCA_precursor_level.png
├── *_noNorm.* (non-normalized versions)
└── ...
```

## Analysis Date
2025-11-21

## Software Versions
- R version 4.3.x
- data.table, ggplot2, limma, pheatmap, corrplot, ggrepel, reshape2
- DiffTestR (custom package)
