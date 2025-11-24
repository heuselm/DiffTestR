#' Import precursor-level matrix from .tsv table (wide format)
#'
#' @description Import precursor-level .tsv table (wide format)
#'
#' @param path_to_tsv Path to wide format result tsv file, with column
#' names Protein.Group, Precursor.Id, filename1 filename2 etc. Can also include
#' annotation columns like Protein.Ids, Protein.Names, Genes, etc.
#' @param precursor_annotation Optional tsv format table with columns Precursor.Id,
#' Protein.Group and extended information. If NULL (default), annotation columns
#' will be extracted from path_to_tsv if available (e.g., Protein.Ids, Protein.Names,
#' Genes, First.Protein.Description, Proteotypic, Stripped.Sequence, Modified.Sequence,
#' Precursor.Charge). This parameter is typically not needed for DIA-NN output files.
#' @param study_design Study design in tab-separated .txt with mandatory
#' columns filename, condition, replicate
#'
#' @return A qdata result object (list) containing objects
#'data, data_wide, matrix_raw, ann_col, ann_row, study_design and source_file
#'
#' @author Moritz Heusel
#'
#' @import data.table
#' @export
#'
importFromTables <- function(path_to_tsv = "../../04_DIANN_LfMS1/Table_precursor_wide.txt",
                              precursor_annotation = NULL,
                              study_design = "../../04_DIANN_LfMS1/Study_design.txt")
{
  # generate container object
  qdata = list()

  # import
  data = data.table::fread(path_to_tsv)
  # data = data.table::fread("../../04_DIANN_LfMS1/DIANN_LfMS1-Pass2.pr_matrix.tsv")
  # names(data)
  # data = data[, c(1,10,11:30)]

  # Clean column names: remove paths and file extensions
  names(data) = sapply(names(data), function(x) {
    # Remove file extensions (.raw, .d, etc.)
    x = gsub("\\.(raw|d|mzML|wiff)$", "", x, ignore.case = TRUE)
    # Extract basename (handle both Windows \ and Unix / paths)
    x = basename(x)
    return(x)
  })

  # keep only data for files in Study design
  study_des = data.table::fread(study_design)
  study_des$filename = gsub("\\.(raw|d|mzML|wiff)$", "",
                            basename(study_des$filename),
                            ignore.case = TRUE)

  data.s = data[, c(which(names(data) %in% "Protein.Group"),
                    which(names(data) %in% "Precursor.Id"),
                    which(names(data) %in% study_des$filename)), with = F]

  if (nrow(data.s) == 0){
    message("filenames in data: ", cat(unique(names(data)[3:ncol(data)], "\n")))
    message("filenames in ", study_design, ": ", cat(unique(study_des$filename)))
    stop("All runs excluded during import. Make sure file names in study_design match those in the data!")
  }

  # Reformat to long and add study design information
  data.s.long = melt(data.s, id.vars = 1:2, measure.vars = 3:ncol(data.s),
                     variable.name = "File.Name",
                     value.name = "Precursor.Quantity")

  data.s.long = merge(data.s.long, study_des, by.x = "File.Name", by.y = "filename")


  # Reformat to wide, save as matrix and assemble matrix row and column annotations
  data.s.wide = data.table::dcast(data.s.long,
                                  Protein.Group+Precursor.Id~File.Name+condition+replicate,
                                  value.var = "Precursor.Quantity",
                                  fun.aggregate = sum)
  data.s.wide.m.raw = as.matrix(data.s.wide[, 3:ncol(data.s.wide)])
  row.names(data.s.wide.m.raw) = data.s.wide$Precursor.Id

  # assemble annotations
  ann_col = study_des
  ann_col[, rn:=paste0(filename,condition,replicate, sep = "_"), filename]
  row.names(ann_col) = ann_col$rn
  ann_col$rn = NULL

  # Handle precursor annotation
  if (!is.null(precursor_annotation)) {
    # Use external annotation file
    prec_ann = fread(precursor_annotation)
    prec_ann = prec_ann[Precursor.Id != ""]
    ann_row = merge(data.s.wide[,1:2], prec_ann,
                    by = "Precursor.Id", all.x = T, all.y = F)
  } else {
    # Extract annotation columns from the main data file
    # Common DIA-NN annotation columns
    ann_cols = c("Protein.Group", "Precursor.Id", "Protein.Ids", "Protein.Names",
                 "Genes", "First.Protein.Description", "Proteotypic",
                 "Stripped.Sequence", "Modified.Sequence", "Precursor.Charge")

    # Keep only columns that exist in the data
    available_ann_cols = intersect(ann_cols, names(data))

    if (length(available_ann_cols) > 0) {
      # Extract annotation from original data
      ann_row = unique(data[, ..available_ann_cols])
      # Merge with data.s.wide to ensure we only keep precursors that passed filtering
      ann_row = merge(data.s.wide[,1:2], ann_row,
                      by = c("Protein.Group", "Precursor.Id"),
                      all.x = T, all.y = F)
    } else {
      # Minimal annotation if no columns found
      ann_row = data.s.wide[,1:2]
    }
  }

  row.names(ann_row) = data.s.wide$Precursor.Id

  # collect data in list result object
  qdata$data_long = data.s.long
  qdata$data_wide = data.s.wide
  qdata$matrix_raw = data.s.wide.m.raw
  qdata$ann_col = ann_col
  qdata$ann_row = ann_row
  qdata$study_design = study_design
  qdata$source_file = path_to_tsv
  return(qdata)
}
