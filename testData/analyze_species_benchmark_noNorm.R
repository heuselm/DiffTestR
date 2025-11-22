# Species Benchmark Analysis Script - WITHOUT NORMALIZATION
# Analyze HYE (Human-Yeast-E.coli) species mix data from Exploris480 and TimsTOF
# Compare MixA vs MixB for both instruments
# This version skips quantile normalization to analyze raw quantitative values

# Load required packages
library(data.table)
library(ggplot2)
library(gridExtra)
library(limma)
library(pheatmap)
library(corrplot)
library(ggrepel)
library(reshape2)

# Source DiffTestR functions
source("../R/testDifferentialAbundance.R")
source("../R/importFromTables.R")

# Set working directory
setwd("c:/Users/heuse/My Drive/Rcode/DiffTestR/testData")

# Create output directory
if(!dir.exists("results")) dir.create("results")

#===============================================================================
# Function to extract species from Protein.Names
#===============================================================================
extract_species <- function(protein_names) {
  ifelse(grepl("_HUMAN", protein_names) & !grepl(";", protein_names), "H. sapiens",
  ifelse(grepl("_YEAST", protein_names) & !grepl(";", protein_names), "S. cerevisiae",
  ifelse(grepl("_ECOLI", protein_names) & !grepl(";", protein_names), "E. coli", "Mixed")))
}

#===============================================================================
# Function to create study design from stats file
#===============================================================================
create_study_design <- function(stats_file, instrument_name) {
  stats <- fread(stats_file)

  # Extract sample names from File.Name
  filenames <- basename(stats$File.Name)
  filenames <- gsub("\\.raw$|\\.d$", "", filenames)

  # Determine condition from filename
  # For Exploris: MixA or MixB
  # For TimsTOF: OT2_A or OT2_B
  if (grepl("Exploris", instrument_name)) {
    condition <- ifelse(grepl("MixA", filenames), "A", "B")
  } else {
    condition <- ifelse(grepl("_A_", filenames), "A", "B")
  }

  # Extract replicate number
  replicate <- as.integer(gsub(".*_TR(\\d+).*|.*_R(\\d+).*", "\\1\\2", filenames))

  study_design <- data.table(
    filename = filenames,
    condition = condition,
    replicate = replicate,
    instrument = instrument_name
  )

  return(study_design)
}

#===============================================================================
# EXPLORIS 480 ANALYSIS - NO NORMALIZATION
#===============================================================================
cat("Processing Exploris480 dataset (no normalization)...\n")

# Create study design for Exploris
exploris_design <- create_study_design("Exploris480.stats.tsv", "Exploris480")
fwrite(exploris_design, "results/Exploris480_study_design_noNorm.tsv", sep = "\t")

# Import and process Exploris data
exploris_data <- fread("Exploris480.tsv")

# Clean filenames to match study design
exploris_data[, Run := gsub("\\.raw$", "", basename(Run))]

# Filter to proteins from single species only (discard proteins with ; in Protein.Group)
exploris_data_clean <- exploris_data[!grepl(";", Protein.Group)]

cat(sprintf("Exploris: Removed %d multi-species protein groups\n",
            nrow(exploris_data) - nrow(exploris_data_clean)))

# Add species annotation
exploris_data_clean[, Species := ifelse(grepl("_HUMAN", Protein.Names), "H. sapiens",
  ifelse(grepl("_YEAST", Protein.Names), "S. cerevisiae",
  ifelse(grepl("_ECOLI", Protein.Names), "E. coli", "Unknown")))]

# Create wide format for testDifferentialAbundance
exploris_wide <- dcast(exploris_data_clean,
                       Protein.Group + Protein.Names + Precursor.Id + Species ~ Run,
                       value.var = "Precursor.Quantity")

# Run differential abundance testing WITHOUT NORMALIZATION
exploris_results <- testDifferentialAbundance(
  input_dt = exploris_wide,
  study_design = exploris_design,
  condition_1 = "A",
  condition_2 = "B",
  min_n_obs = 4,
  normalize_data = FALSE,  # CHANGED: No normalization
  imp_percentile = 0.001,
  imp_sd = 0.2,
  plot_pdf = FALSE,
  write_tsv_tables = FALSE
)

#===============================================================================
# TIMSTOF ANALYSIS - NO NORMALIZATION
#===============================================================================
cat("\nProcessing TimsTof dataset (no normalization)...\n")

# Create study design for TimsTof
timstof_design <- create_study_design("TimsTof.stats.tsv", "TimsTof")
fwrite(timstof_design, "results/TimsTof_study_design_noNorm.tsv", sep = "\t")

# Import and process TimsTof data
timstof_data <- fread("TimsTof.tsv")

# Clean filenames to match study design
timstof_data[, Run := gsub("\\.d$", "", basename(Run))]

# Filter to proteins from single species only
timstof_data_clean <- timstof_data[!grepl(";", Protein.Group)]

cat(sprintf("TimsTof: Removed %d multi-species protein groups\n",
            nrow(timstof_data) - nrow(timstof_data_clean)))

# Add species annotation
timstof_data_clean[, Species := ifelse(grepl("_HUMAN", Protein.Names), "H. sapiens",
  ifelse(grepl("_YEAST", Protein.Names), "S. cerevisiae",
  ifelse(grepl("_ECOLI", Protein.Names), "E. coli", "Unknown")))]

# Create wide format for testDifferentialAbundance
timstof_wide <- dcast(timstof_data_clean,
                      Protein.Group + Protein.Names + Precursor.Id + Species ~ Run,
                      value.var = "Precursor.Quantity")

# Run differential abundance testing WITHOUT NORMALIZATION
timstof_results <- testDifferentialAbundance(
  input_dt = timstof_wide,
  study_design = timstof_design,
  condition_1 = "A",
  condition_2 = "B",
  min_n_obs = 4,
  normalize_data = FALSE,  # CHANGED: No normalization
  imp_percentile = 0.001,
  imp_sd = 0.2,
  plot_pdf = FALSE,
  write_tsv_tables = FALSE
)

#===============================================================================
# PREPARE RESULTS FOR VISUALIZATION
#===============================================================================

# Get protein-level results
exploris_protein <- unique(exploris_results$diffExpr_result_dt[,
  .(Protein.Group, Protein.Names, log2_fold_change_protein,
    p_value_BHadj_protein, n_precursors)])

# Add Species information back
exploris_species <- unique(exploris_data_clean[, .(Protein.Group, Species)])
exploris_protein <- merge(exploris_protein, exploris_species, by = "Protein.Group")
exploris_protein[, Instrument := "Exploris 480"]

timstof_protein <- unique(timstof_results$diffExpr_result_dt[,
  .(Protein.Group, Protein.Names, log2_fold_change_protein,
    p_value_BHadj_protein, n_precursors)])

# Add Species information back
timstof_species <- unique(timstof_data_clean[, .(Protein.Group, Species)])
timstof_protein <- merge(timstof_protein, timstof_species, by = "Protein.Group")
timstof_protein[, Instrument := "timsTOF Pro 2"]

# Calculate mean abundance for MA plot
# Sum precursor intensities per protein
exploris_abundance <- exploris_data_clean[, .(mean_abundance = mean(log10(Precursor.Quantity + 1), na.rm = TRUE)),
                                          by = .(Protein.Group, Species)]
exploris_protein <- merge(exploris_protein, exploris_abundance, by = c("Protein.Group", "Species"))

timstof_abundance <- timstof_data_clean[, .(mean_abundance = mean(log10(Precursor.Quantity + 1), na.rm = TRUE)),
                                        by = .(Protein.Group, Species)]
timstof_protein <- merge(timstof_protein, timstof_abundance, by = c("Protein.Group", "Species"))

# Combine both datasets
combined_results <- rbind(exploris_protein, timstof_protein)

# Add significance category
combined_results[, Significant := p_value_BHadj_protein <= 0.05 & abs(log2_fold_change_protein) >= 1]

# Save combined results with _noNorm suffix
fwrite(combined_results, "results/combined_protein_results_noNorm.tsv", sep = "\t")

#===============================================================================
# VOLCANO PLOTS - Side by side, colored by species
#===============================================================================
cat("\nCreating volcano plots (no normalization)...\n")

# Define colors for species
species_colors <- c("H. sapiens" = "#E41A1C",
                   "S. cerevisiae" = "#377EB8",
                   "E. coli" = "#4DAF4A")

# Create volcano plot
volcano_plot <- ggplot(combined_results,
                      aes(x = log2_fold_change_protein,
                          y = -log10(p_value_BHadj_protein),
                          color = Species)) +
  geom_point(alpha = 0.6, size = 1.5) +
  facet_wrap(~Instrument, ncol = 2) +
  scale_color_manual(values = species_colors) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray40") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray40") +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.3) +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.3) +
  labs(x = expression(paste("Fold change, ", log[2])),
       y = expression(paste("LFQ intensity, ", -log[10], "p")),
       title = "Differential Abundance Analysis: Mix A vs Mix B (No Normalization)",
       subtitle = "Species-specific protein quantification") +
  theme_bw() +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(face = "bold", size = 11),
    panel.grid.minor = element_blank()
  )

ggsave("results/volcano_plots_sidebyside_noNorm.pdf", volcano_plot, width = 12, height = 6)
ggsave("results/volcano_plots_sidebyside_noNorm.png", volcano_plot, width = 12, height = 6, dpi = 300)

#===============================================================================
# MA PLOTS - A/B ratio over summed log10 abundance
#===============================================================================
cat("Creating MA plots (ratio vs abundance, no normalization)...\n")

ma_plot <- ggplot(combined_results,
                 aes(x = mean_abundance,
                     y = log2_fold_change_protein,
                     color = Species)) +
  geom_point(alpha = 0.6, size = 1.5) +
  facet_wrap(~Instrument, ncol = 2) +
  scale_color_manual(values = species_colors) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 0.5) +
  geom_hline(yintercept = c(-1, 1), linetype = "dashed", color = "gray40") +
  labs(x = expression(paste("Summed ", log[10], " abundance")),
       y = expression(paste("A/B ratio, ", log[2])),
       title = "MA Plot: Mix A vs Mix B (No Normalization)",
       subtitle = "Ratio vs abundance colored by species") +
  theme_bw() +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(face = "bold", size = 11),
    panel.grid.minor = element_blank()
  )

ggsave("results/ma_plots_sidebyside_noNorm.pdf", ma_plot, width = 12, height = 6)
ggsave("results/ma_plots_sidebyside_noNorm.png", ma_plot, width = 12, height = 6, dpi = 300)

#===============================================================================
# SUMMARY STATISTICS
#===============================================================================
cat("\n=== SUMMARY STATISTICS (NO NORMALIZATION) ===\n")

summary_stats <- combined_results[, .(
  n_proteins = .N,
  n_significant = sum(Significant),
  median_log2FC = median(log2_fold_change_protein, na.rm = TRUE),
  median_pvalue = median(p_value_BHadj_protein, na.rm = TRUE)
), by = .(Instrument, Species)]

print(summary_stats)
fwrite(summary_stats, "results/summary_statistics_noNorm.tsv", sep = "\t")

# Species-specific fold changes
cat("\n=== Expected vs Observed Fold Changes (NO NORMALIZATION) ===\n")
fc_summary <- combined_results[, .(
  median_log2FC = median(log2_fold_change_protein, na.rm = TRUE),
  mean_log2FC = mean(log2_fold_change_protein, na.rm = TRUE),
  sd_log2FC = sd(log2_fold_change_protein, na.rm = TRUE)
), by = .(Instrument, Species)]

print(fc_summary)

cat("\nAnalysis complete! Results saved in 'results/' directory.\n")
cat("Generated files:\n")
cat("  - volcano_plots_sidebyside_noNorm.pdf\n")
cat("  - ma_plots_sidebyside_noNorm.pdf\n")
cat("  - combined_protein_results_noNorm.tsv\n")
cat("  - summary_statistics_noNorm.tsv\n")
