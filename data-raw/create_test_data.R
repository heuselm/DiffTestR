# Create test dataset from Exploris480 data
# Subset to 300 proteins for package testing

library(data.table)

# Load full Exploris480 dataset
exploris_full <- fread("testData/Exploris480.tsv")

# Set seed for reproducibility
set.seed(42)

# Select 300 random proteins
# Ensure we include proteins from different species
human_proteins <- unique(exploris_full[grepl("_HUMAN", Protein.Names) & !grepl(";", Protein.Group), Protein.Group])
yeast_proteins <- unique(exploris_full[grepl("_YEAST", Protein.Names) & !grepl(";", Protein.Group), Protein.Group])
ecoli_proteins <- unique(exploris_full[grepl("_ECOLI", Protein.Names) & !grepl(";", Protein.Group), Protein.Group])

# Sample 100 proteins from each species
set.seed(42)
selected_proteins <- c(
  sample(human_proteins, min(100, length(human_proteins))),
  sample(yeast_proteins, min(100, length(yeast_proteins))),
  sample(ecoli_proteins, min(100, length(ecoli_proteins)))
)

# Filter data to selected proteins
test_data_long <- exploris_full[Protein.Group %in% selected_proteins]

cat("Created test dataset with:\n")
cat("  - Proteins:", length(unique(test_data_long$Protein.Group)), "\n")
cat("  - Precursors:", length(unique(test_data_long$Precursor.Id)), "\n")
cat("  - Observations:", nrow(test_data_long), "\n")

# Get actual filenames from the data
actual_filenames <- unique(test_data_long$Run)
cat("Actual filenames in data:\n", actual_filenames, "\n\n")

# Create study design matching actual filenames
test_study_design <- data.table(
  filename = actual_filenames,
  condition = ifelse(grepl("MixA", actual_filenames), "A", "B"),
  replicate = as.integer(gsub(".*TR(\\d+).*", "\\1", actual_filenames))
)
test_study_design <- test_study_design[order(condition, replicate)]

# Clean filenames in test data
test_data_long[, Run := gsub("\\.raw$", "", basename(Run))]

# Create wide format for testing
test_data_wide <- dcast(test_data_long,
                        Protein.Group + Protein.Names + Precursor.Id ~ Run,
                        value.var = "Precursor.Quantity")

# Save as R data objects
DiffTestR_example_data_long <- test_data_long
DiffTestR_example_data_wide <- test_data_wide
DiffTestR_example_study_design <- test_study_design

# Save to data/ directory
usethis::use_data(DiffTestR_example_data_long, overwrite = TRUE)
usethis::use_data(DiffTestR_example_data_wide, overwrite = TRUE)
usethis::use_data(DiffTestR_example_study_design, overwrite = TRUE)

cat("\nTest data saved to data/ directory\n")
cat("Access via:\n")
cat("  - data(DiffTestR_example_data_long)\n")
cat("  - data(DiffTestR_example_data_wide)\n")
cat("  - data(DiffTestR_example_study_design)\n")
