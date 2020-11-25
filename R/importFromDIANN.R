#' Import precursor-level matrix from DIANN results
#'
#' @description Import precursor-level matrix from DIANN results
#'
#' @param path_to_result Path to DIANN result tsv file, long format
#' @param quant_value_col The name of the column containing precursor level quant
#' values to use. Typically Precursor.Quantity or Precursor.Normalised
#' @param study_design Study design in tab-separated .txt with mandatory columns filename, condition, replicate
#'
#' @return A qdata result object (list) containing objects data, data_wide, matrix_raw, ann_col and ann_row
#'
#' @author Moritz Heusel
#'
#' @import  data.table
#' @export
#'
importFromDIANN <- function(path_to_result = "../../04_DIANN_LfMS1/DIANN_LfMS1-Pass2.tsv",
                              quant_value_col = "Precursor.Quantity",
                              study_design = "../../04_DIANN_LfMS1/Study_design.txt")
{
  # generate container object
  qdata = list()

  # import
  data = data.table::fread(path_to_result)

  # flatten filenames
  data[, File.Name:=unlist(strsplit(File.Name, split = "\\\\"))[length(unlist(strsplit(File.Name, split = "\\\\")))], File.Name]
  data[, File.Name:=gsub(".raw", "", File.Name)]

  # keep only data for files in Study design
  study_des = data.table::fread(study_design)

  # study_des$filename = paste(study_des$filename, "012")

  data.s = data[File.Name %in% study_des$filename]

  if (nrow(data.s) == 0){
    message("filenames in data: ", cat(unique(data$File.Name)), "\n")
    message("filenames in ", study_design, ": ", cat(unique(study_des$filename)))
    stop("All runs excluded during import. Make sure file names in study_design match those in the data!")
  }

  # Add study design information
  data.s = merge(data.s, study_des, by.x = "File.Name", by.y = "filename")

  # Reformat to wide, save as matrix and assemble matrix row and column annotations
  data.s.wide = data.table::dcast(data.s, Protein.Group+Protein.Ids+Protein.Names+Genes+Modified.Sequence+Stripped.Sequence+Precursor.Id~File.Name+condition+replicate, value.var = quant_value_col)
    data.s.wide.m.raw = as.matrix(data.s.wide[, 8:ncol(data.s.wide)])
  row.names(data.s.wide.m.raw) = data.wide$Precursor.Id

  # assemble annotations
  ann_col = study_design
  ann_col[, rn:=paste0(filename,condition,replicate, sep = "_"), filename]
  row.names(ann_col) = ann_col$rn
  ann_col$rn = NULL

  ann_row = data.wide[,1:6]
  row.names(ann_row) = data.wide$Precursor.Id

  # collect data in list result object
  qdata$data_long = data.s
  qdata$data_wide = data.s.wide
  qdata$matrix_raw = data.s.wide.m.raw
  qdata$ann_col = ann_col
  qdata$ann_row = ann_row
  qdata$study_design = study_design
  qdata$source_file = path_to_result
  return(qdata)
}
