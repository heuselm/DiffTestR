# differential expression testing on precursor level
# M.H. -> Tanner2021
#######################################################

# load packages/dependencies
library(devtools)
devtools::install_github("heuselm/DiffTestR", ref = "Tanner2021")
# for code details,
# see https://github.com/heuselm/DiffTestR/blob/Tanner2021/R/testDifferentialAbundance.R
library(DiffTestR)
library(data.table)
library(ggplot2)
library(ggrepel)

# load input data
input_long = fread("PXD029625/DIANN_results_Lung/210331_R04_Lung_DIANN1712_LibFree2Pass/Pass2.tsv")
input_long[, File.Name:=unlist(strsplit(File.Name, split = "\\\\"))
  [length(unlist(strsplit(File.Name, split = "\\\\")))], File.Name]
input_long[, File.Name:=gsub(".raw", "", File.Name)]
input_wide = dcast(input_long,
  Precursor.Id+Modified.Sequence+Protein.Group+Protein.Names+Genes~File.Name,
                   value.var = "Precursor.Quantity")

study_design = fread("annotation_lung.csv")
study_design

# Conditions in dataset:
conditions = unique(study_design$condition)
conditions

# Conditions of interest for testing & visualization:
pairs = matrix(c("BTH", "DEX",
  "TH", "V",
  "B", "B",
  "B", "B"), ncol =2)
pairs

cRes = data.table()

# run differential expression testing for these comparisons
for (i in 1:nrow(pairs)){
  c1 = pairs[i,1]
  c2 = pairs[i,2]

  dir.create(paste0(i, "_", c1, "_vs_", c2, "_Lung"))
  setwd(paste0(i, "_", c1, "_vs_", c2, "_Lung"))

  diffTestRes = testDifferentialAbundance(input_dt = input_wide, 
                                          study_design = study_design,
                                          condition_1 = c1,
                                          condition_2 = c2,
                                          # filtering for min. number of observations
                                          min_n_obs = 4,
                                          # imputation of missing values
                                          imp_percentile = 0.001,
                                          imp_sd = 0.2,
                                          # plots?
                                          plot_pdf = TRUE,
                                          # tsv result table?
                                          write_tsv_tables = TRUE,
                                          # highlight protein in volcano?
                                          target_protein = "A0A571BF69")

  saveRDS(diffTestRes, "diffTestRes.rds")
  cRes = rbind(cRes, diffTestRes$diffExpr_result_dt[,
    .(Protein.Group, Protein.Names, comparison,
      p_value, log2_fold_change,target_prot,
      log2_fold_change_protein,n_precursors,p_value_protein,p_value_BHadj_protein)])

  setwd("../")
}

names(cRes)

# Visualize combined results
DiffTestR::plotDifferentialAbundanceOverview(cRes, label_prefix = "01_diffExpr_Lung", protein_highlight_tag = NA)

# export combined results
fwrite(cRes[, compartment:="Lung"], "differentialTestingResults_Lung.csv")