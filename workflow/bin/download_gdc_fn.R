#' Download and Save GDC Gene Expression Data
#'
#' This function queries, downloads, and saves gene expression quantification
#' data from the Genomic Data Commons (GDC) for a specified cancer project.
#'
#' @param cancer Character string specifying the GDC project
#'   (e.g., "TCGA-BRCA").
#' @param output_file Character string specifying the file path where the output
#'   data will be saved.
#'
#' @return Saves an RData file containing the gene expression data in the
#'   specified output file.
#' @details
#' The function uses the `GDCquery`, `GDCdownload`, and `GDCprepare` functions
#' from the TCGAbiolinks package to retrieve and process gene expression
#' quantification data (STAR - Counts workflow) for the given cancer project.
#'
#' @importFrom TCGAbiolinks GDCquery
#' @importFrom TCGAbiolinks GDCdownload
#' @importFrom TCGAbiolinks GDCprepare
download_gdc_expression <- function(cancer, output_file, files_per_chunk = 7) {
  mrna_query <- GDCquery(
    project = cancer,
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "STAR - Counts"
  )
  cat("downloading data", "\n")
  GDCdownload(mrna_query, method = "api", files.per.chunk = files_per_chunk)
  cat("downloading data", "\n")
  mrna_df <- GDCprepare(mrna_query)
  cat("saving expression data", "\n")
  save(mrna_df,
       file = output_file)
}


#' Download and Save GDC Clinical Data
#'
#' This function downloads clinical data for a specified cancer project and
#' saves it as an RData file in the given output file.
#'
#' @param cancer Character string specifying the GDC project
#'   (e.g., "TCGA-BRCA").
#' @param output_file Character string specifying the file path where the output
#'   data will be saved.
#'
#' @return Saves an RData file with clinical data in the specified file.
#' @importFrom TCGAbiolinks GDCquery_clinic
download_gdc_clinical <- function(cancer, output_file) {
    clinical <- GDCquery_clinic(
        project = cancer,
        type = "clinical",
        save.csv = FALSE
    )
    save(
        clinical,
        file = output_file
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

#' Read Summarized Experiment Data
#'
#' Loads a SummarizedExperiment object and extracts its assay data.
#'
#' @param data_file Path to the RData file with a SummarizedExperiment 'mrna_df'.
#' @return Matrix or data frame with assay data from the SummarizedExperiment.
#' @import SummarizedExperiment
#' @export
readSummarizedExperiment <- function(data_file) {
    load(data_file, exp <- new.env())
    mrna_df <- exp[['mrna_df']]
    mrna_df <- SummarizedExperiment::assay(mrna_df)
    return(mrna_df)
}