# CV and PCA Analysis for Species Benchmark
# Compare coefficient of variation at precursor vs protein level
# Create PCA plots showing A vs B separation per instrument

# Load required packages
library(data.table)
library(ggplot2)
library(gridExtra)

# Set working directory
setwd("c:/Users/heuse/My Drive/Rcode/DiffTestR/testData")

# Create output directory if it doesn't exist
if(!dir.exists("results")) dir.create("results")

#===============================================================================
# Function to calculate CV (coefficient of variation)
#===============================================================================
calc_cv <- function(x) {
  sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE) * 100
}

#===============================================================================
# EXPLORIS 480 - CV ANALYSIS
#===============================================================================
cat("Processing Exploris480 CV analysis...\n")

# Load data
exploris_data <- fread("Exploris480.tsv")
exploris_data[, Run := gsub("\\.raw$", "", basename(Run))]

# Filter to single species proteins
exploris_data_clean <- exploris_data[!grepl(";", Protein.Group)]

# Add species annotation
exploris_data_clean[, Species := ifelse(grepl("_HUMAN", Protein.Names), "H. sapiens",
  ifelse(grepl("_YEAST", Protein.Names), "S. cerevisiae",
  ifelse(grepl("_ECOLI", Protein.Names), "E. coli", "Unknown")))]

# Add condition labels
exploris_data_clean[, Condition := ifelse(grepl("MixA", Run), "A", "B")]

# Calculate CV at PRECURSOR level (within each condition)
exploris_precursor_cv <- exploris_data_clean[!is.na(Precursor.Quantity), .(
  CV_precursor = calc_cv(Precursor.Quantity),
  mean_abundance = mean(Precursor.Quantity, na.rm = TRUE),
  n_obs = sum(!is.na(Precursor.Quantity))
), by = .(Precursor.Id, Protein.Group, Species, Condition)]

# Filter to precursors with at least 3 observations
exploris_precursor_cv <- exploris_precursor_cv[n_obs >= 3]

# Calculate CV at PROTEIN level (within each condition)
# First aggregate to protein level per run
exploris_protein_quant <- exploris_data_clean[!is.na(Precursor.Quantity), .(
  Protein.Quantity = sum(Precursor.Quantity, na.rm = TRUE)
), by = .(Protein.Group, Species, Run, Condition)]

exploris_protein_cv <- exploris_protein_quant[, .(
  CV_protein = calc_cv(Protein.Quantity),
  mean_abundance = mean(Protein.Quantity, na.rm = TRUE),
  n_obs = .N
), by = .(Protein.Group, Species, Condition)]

# Filter to proteins with at least 3 observations
exploris_protein_cv <- exploris_protein_cv[n_obs >= 3]

# Add instrument label
exploris_precursor_cv[, Instrument := "Exploris 480"]
exploris_protein_cv[, Instrument := "Exploris 480"]

#===============================================================================
# TIMSTOF - CV ANALYSIS
#===============================================================================
cat("Processing TimsTof CV analysis...\n")

# Load data
timstof_data <- fread("TimsTof.tsv")
timstof_data[, Run := gsub("\\.d$", "", basename(Run))]

# Filter to single species proteins
timstof_data_clean <- timstof_data[!grepl(";", Protein.Group)]

# Add species annotation
timstof_data_clean[, Species := ifelse(grepl("_HUMAN", Protein.Names), "H. sapiens",
  ifelse(grepl("_YEAST", Protein.Names), "S. cerevisiae",
  ifelse(grepl("_ECOLI", Protein.Names), "E. coli", "Unknown")))]

# Add condition labels (TimsTOF uses _A_ and _B_ pattern)
timstof_data_clean[, Condition := ifelse(grepl("_A_", Run), "A", "B")]

# Calculate CV at PRECURSOR level (within each condition)
timstof_precursor_cv <- timstof_data_clean[!is.na(Precursor.Quantity), .(
  CV_precursor = calc_cv(Precursor.Quantity),
  mean_abundance = mean(Precursor.Quantity, na.rm = TRUE),
  n_obs = sum(!is.na(Precursor.Quantity))
), by = .(Precursor.Id, Protein.Group, Species, Condition)]

# Filter to precursors with at least 3 observations
timstof_precursor_cv <- timstof_precursor_cv[n_obs >= 3]

# Calculate CV at PROTEIN level (within each condition)
timstof_protein_quant <- timstof_data_clean[!is.na(Precursor.Quantity), .(
  Protein.Quantity = sum(Precursor.Quantity, na.rm = TRUE)
), by = .(Protein.Group, Species, Run, Condition)]

timstof_protein_cv <- timstof_protein_quant[, .(
  CV_protein = calc_cv(Protein.Quantity),
  mean_abundance = mean(Protein.Quantity, na.rm = TRUE),
  n_obs = .N
), by = .(Protein.Group, Species, Condition)]

# Filter to proteins with at least 3 observations
timstof_protein_cv <- timstof_protein_cv[n_obs >= 3]

# Add instrument label
timstof_precursor_cv[, Instrument := "timsTOF Pro 2"]
timstof_protein_cv[, Instrument := "timsTOF Pro 2"]

#===============================================================================
# COMBINE AND COMPARE CV
#===============================================================================
cat("Comparing CV values...\n")

# Combine data
all_precursor_cv <- rbind(exploris_precursor_cv, timstof_precursor_cv)
all_protein_cv <- rbind(exploris_protein_cv, timstof_protein_cv)

# Calculate summary statistics
cv_summary <- rbind(
  all_precursor_cv[, .(
    Level = "Precursor",
    Median_CV = median(CV_precursor, na.rm = TRUE),
    Mean_CV = mean(CV_precursor, na.rm = TRUE),
    SD_CV = sd(CV_precursor, na.rm = TRUE),
    N = .N
  ), by = .(Instrument, Condition)],

  all_protein_cv[, .(
    Level = "Protein",
    Median_CV = median(CV_protein, na.rm = TRUE),
    Mean_CV = mean(CV_protein, na.rm = TRUE),
    SD_CV = sd(CV_protein, na.rm = TRUE),
    N = .N
  ), by = .(Instrument, Condition)]
)

print(cv_summary)
fwrite(cv_summary, "results/CV_summary_statistics.tsv", sep = "\t")

# Create comparison plots
# Violin plots of CV distributions
precursor_cv_plot <- ggplot(all_precursor_cv[CV_precursor <= 100],
                             aes(x = Instrument, y = CV_precursor, fill = Condition)) +
  geom_violin(alpha = 0.7) +
  geom_boxplot(width = 0.2, position = position_dodge(0.9)) +
  labs(title = "Coefficient of Variation at Precursor Level",
       y = "CV (%)", x = "") +
  theme_bw() +
  theme(legend.position = "bottom") +
  facet_wrap(~Species, ncol = 3)

protein_cv_plot <- ggplot(all_protein_cv[CV_protein <= 100],
                           aes(x = Instrument, y = CV_protein, fill = Condition)) +
  geom_violin(alpha = 0.7) +
  geom_boxplot(width = 0.2, position = position_dodge(0.9)) +
  labs(title = "Coefficient of Variation at Protein Level",
       y = "CV (%)", x = "") +
  theme_bw() +
  theme(legend.position = "bottom") +
  facet_wrap(~Species, ncol = 3)

combined_cv_plot <- grid.arrange(precursor_cv_plot, protein_cv_plot, ncol = 1)

ggsave("results/CV_comparison.pdf", combined_cv_plot, width = 12, height = 10)
ggsave("results/CV_comparison.png", combined_cv_plot, width = 12, height = 10, dpi = 300)

#===============================================================================
# PCA ANALYSIS ON PRECURSOR LEVEL
#===============================================================================
cat("Creating PCA plots...\n")

# Prepare Exploris data for PCA
exploris_pca_wide <- dcast(exploris_data_clean[!is.na(Precursor.Quantity)],
                            Precursor.Id ~ Run,
                            value.var = "Precursor.Quantity",
                            fun.aggregate = mean)

# Convert to matrix and log transform
exploris_pca_matrix <- as.matrix(exploris_pca_wide[, -1])
rownames(exploris_pca_matrix) <- exploris_pca_wide$Precursor.Id

# Remove rows with any NA values
exploris_pca_matrix_complete <- exploris_pca_matrix[complete.cases(exploris_pca_matrix), ]
exploris_pca_matrix_log <- log2(exploris_pca_matrix_complete + 1)

# Perform PCA (samples are columns, so transpose)
exploris_pca <- prcomp(t(exploris_pca_matrix_log), scale. = TRUE, center = TRUE)

# Create PCA data frame
exploris_pca_df <- as.data.table(exploris_pca$x)
exploris_pca_df[, Sample := rownames(exploris_pca$x)]
exploris_pca_df[, Condition := ifelse(grepl("MixA", Sample), "Mix A", "Mix B")]
exploris_pca_df[, Instrument := "Exploris 480"]

# Calculate variance explained
exploris_var_explained <- round(100 * exploris_pca$sdev^2 / sum(exploris_pca$sdev^2), 1)

# Prepare TimsTof data for PCA
timstof_pca_wide <- dcast(timstof_data_clean[!is.na(Precursor.Quantity)],
                           Precursor.Id ~ Run,
                           value.var = "Precursor.Quantity",
                           fun.aggregate = mean)

# Convert to matrix and log transform
timstof_pca_matrix <- as.matrix(timstof_pca_wide[, -1])
rownames(timstof_pca_matrix) <- timstof_pca_wide$Precursor.Id

# Remove rows with any NA values
timstof_pca_matrix_complete <- timstof_pca_matrix[complete.cases(timstof_pca_matrix), ]
timstof_pca_matrix_log <- log2(timstof_pca_matrix_complete + 1)

# Perform PCA
timstof_pca <- prcomp(t(timstof_pca_matrix_log), scale. = TRUE, center = TRUE)

# Create PCA data frame
timstof_pca_df <- as.data.table(timstof_pca$x)
timstof_pca_df[, Sample := rownames(timstof_pca$x)]
timstof_pca_df[, Condition := ifelse(grepl("MixA", Sample), "Mix A", "Mix B")]
timstof_pca_df[, Instrument := "timsTOF Pro 2"]

# Calculate variance explained
timstof_var_explained <- round(100 * timstof_pca$sdev^2 / sum(timstof_pca$sdev^2), 1)

# Combine PCA results
combined_pca_df <- rbind(exploris_pca_df, timstof_pca_df, fill = TRUE)

# Create PCA plot
pca_plot <- ggplot(combined_pca_df, aes(x = PC1, y = PC2, color = Condition, label = Sample)) +
  geom_point(size = 4, alpha = 0.8) +
  geom_text(vjust = -1, size = 3, show.legend = FALSE) +
  facet_wrap(~Instrument, scales = "free") +
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(title = "PCA on Precursor Level - Sample Separation by Condition",
       x = sprintf("PC1 (Exploris: %s%%, TimsTOF: %s%%)",
                   exploris_var_explained[1], timstof_var_explained[1]),
       y = sprintf("PC2 (Exploris: %s%%, TimsTOF: %s%%)",
                   exploris_var_explained[2], timstof_var_explained[2]))

# Save separate plots with proper variance labels
exploris_pca_plot <- ggplot(exploris_pca_df, aes(x = PC1, y = PC2, color = Condition, label = Sample)) +
  geom_point(size = 4, alpha = 0.8) +
  geom_text(vjust = -1, size = 3, show.legend = FALSE) +
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(title = "PCA on Precursor Level - Exploris 480",
       x = sprintf("PC1 (%s%% variance)", exploris_var_explained[1]),
       y = sprintf("PC2 (%s%% variance)", exploris_var_explained[2]))

timstof_pca_plot <- ggplot(timstof_pca_df, aes(x = PC1, y = PC2, color = Condition, label = Sample)) +
  geom_point(size = 4, alpha = 0.8) +
  geom_text(vjust = -1, size = 3, show.legend = FALSE) +
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(title = "PCA on Precursor Level - timsTOF Pro 2",
       x = sprintf("PC1 (%s%% variance)", timstof_var_explained[1]),
       y = sprintf("PC2 (%s%% variance)", timstof_var_explained[2]))

pca_combined <- grid.arrange(exploris_pca_plot, timstof_pca_plot, ncol = 2)

ggsave("results/PCA_precursor_level.pdf", pca_combined, width = 14, height = 6)
ggsave("results/PCA_precursor_level.png", pca_combined, width = 14, height = 6, dpi = 300)

#===============================================================================
# SUMMARY AND CONCLUSIONS
#===============================================================================
cat("\n=== CV ANALYSIS SUMMARY ===\n")
cat("\nMedian CV values by instrument and level:\n")
print(cv_summary[, .(Instrument, Condition, Level, Median_CV, Mean_CV)])

# Calculate overall CV by instrument
overall_cv <- rbind(
  all_precursor_cv[, .(
    Level = "Precursor",
    Median_CV = median(CV_precursor, na.rm = TRUE)
  ), by = Instrument],

  all_protein_cv[, .(
    Level = "Protein",
    Median_CV = median(CV_protein, na.rm = TRUE)
  ), by = Instrument]
)

cat("\n=== OVERALL COMPARISON ===\n")
print(overall_cv)

cat("\n=== CONCLUSIONS ===\n")
cat("1. Coefficient of Variation Comparison:\n")
if (overall_cv[Instrument == "timsTOF Pro 2" & Level == "Precursor", Median_CV] >
    overall_cv[Instrument == "Exploris 480" & Level == "Precursor", Median_CV]) {
  cat("   - TimsTOF shows HIGHER variability at precursor level\n")
  cat("   - This increased variability may explain fewer significantly regulated proteins\n")
} else {
  cat("   - TimsTOF shows SIMILAR or LOWER variability at precursor level\n")
  cat("   - Variability alone does not explain lack of regulated proteins\n")
}

cat("\n2. PCA Analysis:\n")
cat(sprintf("   - Exploris 480: PC1 explains %s%% variance\n", exploris_var_explained[1]))
cat(sprintf("   - timsTOF Pro 2: PC1 explains %s%% variance\n", timstof_var_explained[1]))
if (exploris_var_explained[1] > timstof_var_explained[1]) {
  cat("   - Exploris shows better separation between Mix A and Mix B\n")
} else {
  cat("   - Both instruments show comparable separation\n")
}

cat("\nAnalysis complete! Results saved in 'results/' directory.\n")
cat("Generated files:\n")
cat("  - CV_summary_statistics.tsv\n")
cat("  - CV_comparison.pdf/png\n")
cat("  - PCA_precursor_level.pdf/png\n")
