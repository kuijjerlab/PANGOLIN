#' Download and Save GDC Gene Expression Data
#'
#' This function queries, downloads, and saves gene expression quantification
#' data from the Genomic Data Commons (GDC) for a specified cancer project.
#'
#' @param cancer Character string specifying the GDC project
#'   (e.g., "TCGA-BRCA").
#' @param output_dir Character string specifying the directory where the output
#'   file will be saved.
#'
#' @return Saves an RData file containing the gene expression data in the
#'   specified output directory.
#' @details
#' The function uses the `GDCquery`, `GDCdownload`, and `GDCprepare` functions
#' from the TCGAbiolinks package to retrieve and process gene expression
#' quantification data (STAR - Counts workflow) for the given cancer project.
#' 
#' @import TCGAbiolinks
download_gdc_expression <- function(cancer, output_dir) {
mrna_query <- GDCquery(project = cancer,
                      data.category = "Transcriptome Profiling",
                      data.type = "Gene Expression Quantification",
                      workflow.type = "STAR - Counts")
  cat("downloading data", "\n")
  GDCdownload(mrna_query, method = "api")
  cat("downloading data")
  mrna_df <- GDCprepare(mrna_query)
  cat("saving expression data", "\n")
  save(mrna_df,
      file = file.path(output_dir, paste0(cancer, ".RData")))
}


#' Download and Save GDC Clinical Data
#'
#' This function downloads clinical data for a specified cancer project and
#' saves it as an RData file in the given output directory.
#'
#' @param cancer Character string specifying the GDC project 
#'   (e.g., "TCGA-BRCA").
#' @param output_dir Character string specifying the directory to save the file.
#'
#' @return Saves an RData file with clinical data in the specified directory.
#' @import TCGAbiolinks
download_gdc_clinical <- function(cancer, output_dir) {
    clinical <- GDCquery_clinic(
        cancer, 
        type = "clinical", 
        save.csv = FALSE
    )
    save(
        clinical,
        file = file.path(output_dir, paste0(cancer, "_clinical.RData"))
    )
}

# download specimen data

#' Download Biospecimen Data from GDC
#'
#' Query and download biospecimen data for a specified cancer project from GDC.
#'
#' @param cancer Character string specifying the GDC project ID 
#'   (e.g., "TCGA-BRCA").
#' @return No return value. Downloads biospecimen data to the default GDC 
#'   directory.
#' @import TCGAbiolinks

download_gdc_biospecimen <- function(cancer) {
    query_biospecimen <- GDCquery(
        project = cancer,
        data.category = "Biospecimen",
        data.type = "Biospecimen Supplement",
        data.format = "BCR Biotab"
    )
    GDCdownload(query_biospecimen)
}

