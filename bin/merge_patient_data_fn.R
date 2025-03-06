#' Load and Combine PD-1 Scores with PD-L1 Expression
#'
#' Reads PD-1 scores and PD-L1 expression from a directory and merges them.
#'
#' @param tumor_pd1_dir Character. Directory containing PD-1/PD-L1 data.
#' 
#' @return A data frame with PD-1 scores and PD-L1 (CD274) expression.
#' @export
#'
#' @examples
#' pd1_data <- combine_pd1_scores_pdl1_expression("path/to/data")

combine_pd1_scores_pdl1_expression <- function(tumor_pd1_dir) {
        if (!dir.exists(tumor_pd1_dir)) stop("Error: PD-1 directory not found.")

        # Load PD-1 scores
        pd1_scores <- tryCatch({
            load_pd1_generic(tumor_pd1_dir, type = "pd1_scores")
        }, error = function(e) stop("Error loading PD-1 scores: ", e$message))

        pd1_scores <- as.data.frame(t(pd1_scores))
        if (nrow(pd1_scores) == 0) stop("Error: PD-1 scores data is empty.")
        # Load PD-L1 expression
        pdl1_expression <- tryCatch({
            load_pd1_generic(tumor_pd1_dir, type = "pdl1_expression")
        }, error = function(e) 
            stop("Error loading PD-L1 expression: ", e$message))

        if (!"CD274" %in% colnames(pdl1_expression)) {
            stop("Error: Missing 'CD274' column in PD-L1 data.")
        }

        # Merge data by barcode
        pd1_scores$CD274_exp <- pdl1_expression$CD274[
            match(rownames(pd1_scores), pdl1_expression$bcr_patient_barcode)]
        return(pd1_scores)
        }


#' Read Cox Model Predicted Risk Scores
#'
#' Loads risk score data from a file.
#'
#' @param risk_score_file Character. Path to the risk scores file.
#' 
#' @return A data frame with predicted risk scores.
#' @export


read_cox_predicted_risk_scores <- function(risk_score_file) {
        if (!file.exists(risk_score_file)) {
            stop("Error: Risk score file does not exist.")
        }

        risk_scores <- tryCatch({
            fread(risk_score_file)
        }, error = function(e) stop("Error reading risk scores: ", e$message))

        if (nrow(risk_scores) == 0) stop("Error: Risk scores file is empty.")

        return(risk_scores)
        }




#' Preprocess Immune Data
#'
#' Reads and cleans immune data, adding `bcr_patient_barcode` and removing 
#' unnecessary columns.
#'
#' @param immune_file Character. Path to the immune data file.
#' 
#' @return A cleaned data frame with immune data.
#' @export

preprocess_immune <- function(immune_file) {
        if (!file.exists(immune_file)) {
            stop("Error: Immune file does not exist.")
        }

        immune_data <- tryCatch({
            fread(immune_file)
        }, error = function(e) stop("Error reading immune data: ", e$message))

        if (nrow(immune_data) == 0) stop("Error: Immune data file is empty.")

        # Generate patient barcode
        immune_data$bcr_patient_barcode <- make_bcr_code(immune_data$Mixture)

        # Remove unnecessary columns
        cols_to_remove <- c("P-value", "Correlation", "RMSE", "Mixture")
        cols_to_remove <- intersect(cols_to_remove, colnames(immune_data))
        immune_data <- immune_data[, !colnames(immune_data) %in% cols_to_remove, with = FALSE]
        return(immune_data)
    }

#' Merge Patient Data with PD-1, Risk Scores, and Immune Data
#'
#' This function reads and merges PD-1 expression, risk scores, and immune data, 
#' ensuring all datasets align by `bcr_patient_barcode`.
#'
#' @param tumor_pd1_dir Character. Directory containing tumor PD-1 data.
#' @param risk_score_file Character. Path to the predicted risk scores file.
#' @param immune_file Character. Path to the immune data file.
#' 
#' @return A merged data frame with PD-1 expression, risk scores, and immune data.
#' @export
#' 
merge_patient_data <- function(tumor_pd1_dir, 
                    risk_score_file,
                    immune_file) {
        # Error handling: Check if files exist
        if (!dir.exists(tumor_pd1_dir)) 
            stop("Error: Tumor PD-1 directory does not exist.")
        if (!file.exists(risk_score_file)) 
            stop("Error: Risk score file does not exist.")
        if (!file.exists(immune_file)) 
            stop("Error: Immune file does not exist.")
        # Load and process PD-1 data
        data <- tryCatch({
            combine_pd1_scores_pdl1_expression(tumor_pd1_dir)
        }, error = function(e) 
            stop("Error loading PD-1 expression data: ", e$message))

        # Ensure data has required structure
        if (nrow(data) == 0) stop("Error: PD-1 data is empty.")
        data$bcr_patient_barcode <- rownames(data)
        # Load and process risk score data
        risk_scores <- tryCatch({
            read_cox_predicted_risk_scores(risk_score_file)
        }, error = function(e) 
            stop("Error loading risk score data: ", e$message))
        # Reshape risk scores
        risk_scores <- tryCatch({
            dcast(risk_scores, bcr_patient_barcode ~ component + type,
                value.var = "predicted_risk")
        }, error = function(e)
            stop("Error transforming risk score data: ", e$message))

        # Merge PD-1 and risk scores
        data <- tryCatch({
            merge(data, risk_scores, by = "bcr_patient_barcode")
        }, error = function(e) 
            stop("Error merging PD-1 and risk scores: ", e$message))
        # Load and process immune data
        immune <- tryCatch({
            preprocess_immune(immune_file)
        }, error = function(e) 
            stop("Error loading immune data: ", e$message))

        if (nrow(immune) == 0) 
            stop("Error: Immune data is empty.")
        # Filter immune data to only include matching patients
        immune_clean <- 
            immune %>% filter(bcr_patient_barcode %in% data$bcr_patient_barcode)
        # Merge all data
        data_all <- tryCatch({
            merge(data, immune_clean, by = "bcr_patient_barcode")
        }, error = function(e) 
            stop("Error merging final dataset: ", e$message))

        return(data_all)
}

#' Read and filter Cox proportional hazards model results
#'
#' This function reads a Cox model results file, filters based on a p-value 
#' threshold, and applies cancer-type-specific filtering criteria.
#'
#' @param cox_results_file A string specifying the path to the Cox results file.
#' @param pval_threshold A numeric value for filtering significant results 
#'   (default = 0.05).
#'
#' @return A filtered data frame with an added `component_type` column.
#' @export
#'
#' @import data.table
#' @importFrom dplyr filter mutate
#'
#' @examples
#' \dontrun{
#'   results <- read_all_coxph_results("cox_results.csv", 0.01)
#' }
read_all_coxph_results <- function(cox_results_file,
                        pval_threshold = 0.05) {
        # Check if the file exists
        if (!file.exists(cox_results_file)) {
            stop("Error: Cox results file does not exist.")
        }

        # Load data
        coxph_results <- data.table::fread(cox_results_file)

        # Ensure required columns exist
        required_cols <- c("pvalue", "type", "cancer", "feature")
        missing_cols <- setdiff(required_cols, colnames(coxph_results))
        if (length(missing_cols) > 0) {
            stop("Error: Missing columns -", 
                    paste(missing_cols, collapse = ", "))
        }
        # Filter based on p-value threshold
        coxph_results <- dplyr::filter(coxph_results, pvalue < pval_threshold)

        # Define PFI-specific cancer types
        pfi_cancer <- c("BRCA", "LGG", "PRAD", "READ", "TGCT", "THCA", "THYM")

        # Apply cancer-type filtering
        coxph_results <- coxph_results %>%
            dplyr::filter((type == "PFI" & cancer %in% pfi_cancer) | 
                        (type == "OS" & !cancer %in% pfi_cancer))

        # Create `component_type` column
        coxph_results <- dplyr::mutate(coxph_results, 
                        component_type = paste(feature, type, sep = "_"))

    return(coxph_results)
}

#' Merge patient data from all cancer subdirectories
#'
#' This function scans a directory for patient data files, reads them, and adds 
#' a `cancer` column based on the filename. The data is then merged into a 
#' single data frame.
#'
#' @param cancer_dir A string specifying the directory with patient data files.
#'
#' @return A data frame with merged patient data from all available cancer types.
#' @export
#' 
merge_patient_data_all_cancers <- function(cancer_dir) {
    # Validate input directory
    if (!dir.exists(cancer_dir)) {
        stop("Error: The specified directory does not exist.")
    }

    # List all matching files
    files <- list.files(
        cancer_dir, 
        pattern = "combined_patient_data",
        recursive = TRUE, full.names = TRUE
    )

    if (length(files) == 0) {
        stop("Error: No matching patient data files found in the directory.")
    }

    # Read and merge all data files
    patient_data <- purrr::map_dfr(files, function(file) {
        cancer_type <- 
            stringr::str_extract(basename(file), "(?<=combined_patient_data_)[A-Z]+")
        df <- data.table::fread(file)
        df$cancer <- cancer_type  # Add cancer type column
        return(df)
    })
    return(patient_data)
}
