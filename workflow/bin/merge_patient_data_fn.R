
#' Combine PD-1 Scores with PD-L1 Expression Data
#'
#' Reads PD-1 scores and PD-L1 expression data from separate files and 
#' combines them by matching patient barcodes. The CD274 expression values 
#' are added to the PD-1 scores data frame.
#'
#' @param pd1_scores_file Character. Path to the file containing PD-1 scores.
#' @param pdl1_expression_file Character. Path to the file containing PD-L1 
#'   expression data.
#'
#' @return A data frame with PD-1 scores and corresponding CD274 expression 
#'   values.
#'
#' @export
combine_pd1_scores_pdl1_expression <- function(pd1_scores_file,
                                               pdl1_expression_file) {
        if (!file.exists(pd1_scores_file)) 
            stop("Error: PD-1 scores file not found.")
        if (!file.exists(pdl1_expression_file)) 
            stop("Error: PD-L1 expression file not found.")

        # Load PD-1 scores
        pd1_scores <- tryCatch({
            load_process_pd1_data(pd1_scores_file, type = "pd1_scores")
        }, error = function(e) stop("Error loading PD-1 scores: ", e$message))

        pd1_scores <- as.data.frame(t(pd1_scores))
        if (nrow(pd1_scores) == 0) stop("Error: PD-1 scores data is empty.")
        # Load PD-L1 expression
        pdl1_expression <- tryCatch({
            load_process_pd1_data(pdl1_expression_file,
                                  type = "pdl1_expression")
        }, error = function(e) 
            stop("Error loading PD-L1 expression: ", e$message))

        if (!"CD274" %in% colnames(pdl1_expression)) {
            stop("Error: Missing 'CD274' column in PD-L1 data.")
        }

        # Merge data by barcode
        pd1_scores$CD274_exp <- pdl1_expression$CD274[
            match(rownames(pd1_scores), 
                  pdl1_expression$bcr_patient_barcode)]
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
        immune_data <- immune_data[, !colnames(immune_data) %in% 
                                   cols_to_remove, with = FALSE]
        return(immune_data)
    }

#' Merge Patient Data with PD-1, Risk Scores, and Immune Data
#'
#' This function reads and merges PD-1 expression, risk scores, and immune 
#' data, ensuring all datasets align by `bcr_patient_barcode`.
#'
#' @param pd1_scores_file Character. Path to the PD-1 scores file.
#' @param pdl1_expression_file Character. Path to the PD-L1 expression file.
#' @param risk_score_file Character. Path to the predicted risk scores file.
#' @param immune_file Character. Path to the immune data file.
#' 
#' @return A merged data frame with PD-1 expression, risk scores, and immune 
#'   data.
#' @export
#' 
merge_patient_data <- function(pd1_scores_file,
                    pdl1_expression_file,
                    risk_score_file,
                    immune_file) {
        # Error handling: Check if files exist
        if (!file.exists(pd1_scores_file)) 
            stop("Error: PD-1 scores file does not exist.")
        if (!file.exists(pdl1_expression_file)) 
            stop("Error: PD-L1 expression file does not exist.")
        if (!file.exists(risk_score_file)) 
            stop("Error: Risk score file does not exist.")
        if (!file.exists(immune_file)) 
            stop("Error: Immune file does not exist.")
        # Load and process PD-1 data
        data <- tryCatch({
            combine_pd1_scores_pdl1_expression(pd1_scores_file, 
                                               pdl1_expression_file)
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

#' Merge Patient Data from All Cancer Types
#'
#' Reads and combines patient data files from multiple cancer types into a 
#' single data frame. Each file is expected to contain combined patient data 
#' for a specific cancer type.
#'
#' @param combined_patient_data_files Character vector. Paths to combined 
#'   patient data files for different cancer types.
#'
#' @return A data frame with merged patient data from all cancer types, 
#'   including a 'cancer' column indicating the cancer type.
#'
#' @import data.table
#' @import purrr
#' @import stringr
#' @export
merge_patient_data_all_cancers <- function(combined_patient_data_files) {
    # Validate input files
    if (length(combined_patient_data_files) == 0) {
        stop("Error: No combined patient data files provided.")
    }

    # Read and merge all data files
    patient_data <- 
        purrr::map_dfr(combined_patient_data_files, function(file) {
        cancer_type <- stringr::str_extract(
            basename(file), "(?<=combined_patient_data_)[A-Z]+")
        df <- data.table::fread(file)
        df$cancer <- cancer_type  # Add cancer type column
        return(df)
    })
    return(patient_data)
}
#' Generate Immune Infiltration Table
#'
#' This function computes correlations between principal components
#' and immune infiltration levels across multiple cancer types.
#'
#' @param cox_results_file Path to the Cox proportional hazards results file.
#' @param combined_patient_data_files Character vector. Paths to combined 
#'   patient data files for different cancer types.
#' @param pval_threshold Numeric. Significance threshold for filtering Cox 
#'   results (default: 0.05).
#' @param correlation_method Character. Correlation method to use 
#'   (default: "spearman").
#'
#' @return A data frame containing correlation results between immune cell 
#'   infiltration and PC scores across different cancer types.
#'
#' @import dplyr
#' @import reshape2
#' @import plyr
#' @export

generate_PC_immune_correlation_table <- function(cox_results_file,
                                combined_patient_data_files,
                                pval_threshold = 0.05,
                                correlation_method = "spearman") {
        cox_res <- read_all_coxph_results(cox_results_file,
                                    pval_threshold = pval_threshold)
        combined_data <- 
            merge_patient_data_all_cancers(combined_patient_data_files) 
        cols_cells <- colnames(combined_data)[9:30]
        cor_res_all <- lapply(1:nrow(cox_res), function(i) {
            tumor <- cox_res$cancer[i]
            pc_component <- cox_res$component[i]
            # Select relevant columns
            immune_cols <- c("bcr_patient_barcode", cols_cells)
            pd1_cols <- c("bcr_patient_barcode", pc_component)

            # Filter once and create two separate datasets
            filtered_data <- combined_data %>% filter(cancer == tumor)
            immune_data <- filtered_data %>% select(all_of(immune_cols))
            pd1_scores <- filtered_data %>% select(all_of(pd1_cols))

            # Reshape data
            data_long <- melt(immune_data)
            data_long <- left_join(data_long, pd1_scores,
                            by = "bcr_patient_barcode")
            cor_data <- ddply(data_long, .(variable),
                function(x) cor(x[[pc_component]], x$value, 
                        method = correlation_method, 
                        use = "complete.obs")) 
                # TODO: add p-value and correct for multiple testing
            colnames(cor_data)[2] <- c("corr")
            # Add metadata
            # cor_data <- cor_data %>%
            #     mutate(tumor = tumor, pc_component = pc_component)
            cor_data$tumor_component <- paste(tumor, 
                                        pc_component, 
                                        sep = "_"
                    )
            return(cor_data)
        })
        cor_res_all <- do.call("rbind", cor_res_all)
        cor_res_all$corr <- abs(cor_res_all$corr)
        cor_res_all  <- cor_res_all[!is.na(cor_res_all$corr),]
        return(cor_res_all)
}

#' Clean Cox Results Using PORCUPINE Results
#'
#' Filters Cox proportional hazards results to include only cancer types 
#' that have significant PD-1 signaling pathway activity according to 
#' PORCUPINE analysis results.
#'
#' @param cox_results_file Character. Path to the Cox results file.
#' @param porcupine_results_file Character. Path to the PORCUPINE results file.
#' @param pval_threshold Numeric. P-value threshold for filtering Cox results 
#'   (default: 0.05).
#'
#' @return A filtered data frame containing Cox results for cancer types with 
#'   significant PD-1 pathway activity.
#'
#' @export
clean_cox_results <- function(cox_results_file, 
                              porcupine_results_file,
                              pval_threshold = 0.05) {
    cox_res <- read_all_coxph_results(cox_results_file,
                                pval_threshold = pval_threshold)
    porcupine_res <- filter_porcupine_results_PD1_pathway(
        porcupine_results_file)
    # Filter cox results based on porcupine results
    cox_res_filtered <- cox_res %>%
        filter(cancer %in% porcupine_res$cancer)

    return(cox_res_filtered)
}

#' Filter Porcupine Results for PD-1 Signaling Pathway
#'
#' Reads a Porcupine results file and filters for the 
#' REACTOME_PD_1_SIGNALING pathway.
#'
#' @param porcupine_results_file Path to the Porcupine results file.
#' @return A data.table with entries for the PD-1 signaling pathway.
#' @import data.table
#' @export
#' 
filter_porcupine_results_PD1_pathway <- function(porcupine_results_file) {
    # Load Porcupine results
    porcupine_res <- data.table::fread(porcupine_results_file) 
    # Filter rows matching the PD-1 signaling pathway
    filtered_results <- porcupine_res[
        grep("REACTOME_PD_1_SIGNALING", pathway, fixed = TRUE)
    ]
    
    # Return filtered results
    return(filtered_results)
}