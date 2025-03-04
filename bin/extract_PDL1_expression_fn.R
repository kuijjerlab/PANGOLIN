# library(data.table)
# library(dplyr)
# library(tidyr)

#' Extract Gene Expression for a Specific Gene
#'
#' @param exp_file Path to the gene expression file (CSV/TSV format).
#' @param samples_file Path to the samples file containing sample IDs.
#' @param gene_id The ID of the gene to extract (default: "CD274").
#' @return A data frame with the extracted gene expression values.
#' @export
extract_gene_expression <- function(exp_file, samples_file, gene_id = "CD274") {
        # Check file existence
        if (!file.exists(exp_file))
                stop("The specified `exp_file` does not exist.")
        if (!file.exists(samples_file)) 
                stop("The specified `samples_file` does not exist.")

        # Validate gene_id
        if (!is.character(gene_id) || length(gene_id) != 1)
                stop("`gene_id` must be a single string.")

        # Read expression data
        exp <- fread(exp_file)
        genes <- exp$V1
        exp <- exp[, -1]

        # Read sample metadata
        samples <- fread(samples_file)

        # Assign column names based on sample IDs
        colnames(exp) <- samples$sample_id

        # Find gene index
        idx_gene <- grep(paste0("^", gene_id, "$"), genes)
        if (length(idx_gene) == 0) 
                stop(paste("Gene ID", gene_id, 
                "not found in the expression file."))

        # Extract expression values for the gene
        gene_exp <- as.data.frame(t(exp[idx_gene, ]))
        colnames(gene_exp) <- gene_id
        gene_exp$cancer <- 
                samples$cancer[match(rownames(gene_exp), samples$sample_id)]
        gene_exp$bcr_patient_barcode <- make_bcr_code(rownames(gene_exp))
        return(gene_exp)
        }

#' Filter Gene Expression Data for a Specific Tumor Type
#'
#' @param gene_exp Data frame from `extract_gene_expression()`.
#' @param tumor Tumor type to filter.
#' @return Filtered gene expression data.
#' 
filter_tumor_expression <- function(gene_exp, tumor) {
    if (!"cancer" %in% colnames(gene_exp))
        stop("The `cancer` column is missing in gene_exp.")

    # Standardize cancer labels if needed
    gene_exp$cancer <- gsub("TCGA-", "", gene_exp$cancer)

    # Filter for tumor type
    gene_exp <- gene_exp[gene_exp$cancer == tumor, ]
    return(gene_exp)
}

#' Process Gene Expression Data with Optional Tumor Filtering and Reshaping
#'
#' @param exp_file Path to the gene expression file.
#' @param samples_file Path to the sample metadata file.
#' @param gene_id The gene ID to extract.
#' @param tumor (Optional) Tumor type to filter. If NULL, returns all data.
#' @return A data frame of processed gene expression data.
#' 
process_gene_expression <- function(exp_file, 
                                samples_file,
                                gene_id = "CD274",
                                tumor) {
        if (!is.null(tumor) && !is.character(tumor)) {
                stop("The `tumor` should be provided.")
        }
        gene_exp <- extract_gene_expression(exp_file, samples_file, gene_id)
        gene_exp <- filter_tumor_expression(gene_exp, tumor)
        wide_data <- gene_exp %>%
            pivot_wider(names_from = bcr_patient_barcode, 
                values_from = all_of(gene_id)) %>%
            mutate(gene = gene_id) %>%
            select(gene, everything())
        rownames(wide_data) <- gene_id
        wide_data <- as.data.frame(wide_data[, -c(1:2)])
        return(gene_exp)
}
