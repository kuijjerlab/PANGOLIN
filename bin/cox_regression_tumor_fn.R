#########################################################################
## Functions performing survival analysis on PD1 network edges, PCs 
## from PD1 subnetwork, and PC from PDL1 subnetwork (CD274)
##########################################################################

libraries <- c("data.table", "dplyr", "doParallel",
            "survminer", "survival", "magrittr", "gridExtra",
            "purrr")
lapply(libraries, library, character.only = TRUE)


#' Generate bcr code from Sample ID
#'
#' This function generates a BCR code from a given sample ID.
#' The BCR code is created by concatenating the first three elements of the
#' sample ID, which are separated by hyphens.
#'
#' @param sample_id A character vector of sample ID for example 
#'                  (patient_id (f.ex TCGA-BL-A13I-01A-11R-A277-07))
#' @return A character vector containing the BCR code (f.ex TCGA-BL-A13I)
#' @export
#' 
make_bcr_code <- function(sample_id) {
        if (!is.character(sample_id)) {
            stop("sample_id must be a character vector.")
        }
        # Process sample_id and extract the first three components
        bcr_code <- sapply(strsplit(sample_id, "-"), function(split_id) {
            if (length(split_id) < 3) {
                return(NA)  # Return NA for malformed input
            }
            paste(split_id[1], split_id[2], split_id[3], sep = "-")
        })
        return(bcr_code)
    }


load_clin_curated_tumor <- function(tumor_clin_file_path) {
        # Validate input: Check if the file exists
        if (!file.exists(tumor_clin_file_path)) {
            stop("The specified file path does not exist: ", tumor_file_path)
        }
        clinical <- tryCatch({data.table::fread(tumor_clin_file_path)
        }, error = function(e) {
            stop("Error reading the file: ", e$message)
        })
        return(clinical)
}

type = "pd1_links"
tumor_pd1_dir = "/storage/kuijjerarea/tatiana/PANGOLIN/data_individual_cancers/ACC/pd1_data"
tumor_clin_file <- "/storage/kuijjerarea/tatiana/PANGOLIN/data_individual_cancers/ACC/clinical/curated_clinical_ACC.txt"

load_pd1_generic <- function(tumor_pd1_dir, type) {
        validate_inputs(tumor_pd1_dir, type)
        file_pattern <- determine_pattern(type)
        tumor_file <- find_tumor_file(tumor_pd1_dir, file_pattern)
        data <- load_process_pd1_data(tumor_pd1_dir, tumor_file, type)
        return(data)
}


validate_inputs <- function(tumor_pd1_dir, type) {
        if (!dir.exists(tumor_pd1_dir)) {
                stop("`pd1_dir` doesn't exist.")
        }
        if (!type %in% c("pd1_links", "pd1_net", 
                         "pd1_scores", "cd274_scores",
                         "cd274_expression")) {
                stop("Invalid `type`. Choose valid type.")
        }
}

determine_pattern <- function(type) {
        patterns <- list(pd1_links = "pd1_edges", pd1_net = "pd1_net", 
                         pd1_scores = "pd1_individual_scores", 
                         cd274_scores = "CD274_individual_scores",
                         cd274_expression = "pdl1_expression")
        return(patterns[[type]])
}


find_tumor_file <- function(tumor_pd1_dir, pattern) {
        pattern_files <- list.files(tumor_pd1_dir, pattern = pattern, 
                                full.names = FALSE)
        if (length(pattern_files) == 0) {
                stop("No files found in: ", pd1_dir)
        }

        if (length(pattern_files) == 0) {
                stop("No file found in: ", tumor_pd1_dir)
        }
        if (length(pattern_files) > 1) {
                stop("Several files found: ", tumor_pd1_dir)
        }
        return(pattern_files[1])
}

load_process_pd1_data <- function(tumor_pd1_dir, tumor_file, type) {
        file_path <- file.path(tumor_pd1_dir, tumor_file)
        if (type == "pd1_links") {
                data <- load_pd1_object(file_path, object_name = "links")
                data <- paste(data$reg, data$tar, sep = "_")
        } else if (type == "pd1_net") {
                data <- load_pd1_object(file_path, object_name = "pd1_net")
                colnames(data) <- make_bcr_code(colnames(data))
        } else if (type == "pd1_scores") {
                scores <- load_pd1_object(file_path, object_name = "ind_scores")
                data <- t(scores)
                colnames(data) <- make_bcr_code(colnames(data))
                data <- data[1:2, ] # only take 2 PCs
        } else if (type == "cd274_scores") {
                scores <- load_pd1_object(file_path, object_name = "ind_scores")
                data <- t(scores)
                colnames(data) <- make_bcr_code(colnames(data))
        } else if (type == "cd274_expression") {
                data <- load_pd1_object(file_path, object_name = "pdl1_expression")
                data <- data.frame(data)
        }
        return(data)
}

load_pd1_object <- function(file, object_name) {
        if (!file.exists(file)) {
                stop("The specified file does not exist: ", file)
        }
        # Load the .RData file into a new environment
        env <- new.env()
        load(file, envir = env)
        # Check if the specified object exists in the environment
        if (!exists(object_name, envir = env)) {
                stop("The object '", object_name, 
                     "' was not found in the file.")
        }
        # Retrieve and return the object
        result <- env[[object_name]]
        return(result)
}

#' Combine clinical information from 
#' https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6066282/ and the clinical 
#' information available with TCGA biolinks package.
#' 
#' This function combines clinical information
#' 
#' @param tumor The tumor type 
#' @param clin Clinical data
#' @return Combined clinical data with TCGA subtype information
#' @export
#' 
combine_clins <- function(tumor_clin_file_path) {
        if (!file.exists(tumor_clin_file_path)) {
                stop("The specified file does not exist: ", tumor_clin_file_path)
        }
        clin_tumor <- load_clin_curated_tumor(tumor_clin_file_path)
        if (!"bcr_patient_barcode" %in% colnames(clin_tumor)) {
            stop("The `clin` data frame must contain a column 
                named `bcr_patient_barcode`.")
        }
        # those tumors are not included in TCGAquery_subtype
        tumor <- unique(clin_tumor$type)
        if (!tumor %in% c("DLBC", "LAML", "MESO", "OV", "TGCT", "THYM")) {
            dataSubt <- tryCatch({
                TCGAbiolinks::TCGAquery_subtype(tumor = tumor)
            }, error = function(e) {
                stop("Error querying TCGA subtype data for tumor type: ", 
                    tumor, ". Ensure the tumor type is supported.")
            })
            colnames(dataSubt) <- gsub(" ", "_", colnames(dataSubt))
            dataSubt$patient_barcode <- dataSubt[,1, drop = TRUE]
            dataSubt$patient_barcode <- as.character(dataSubt$patient_barcode)
            dataSubt$bcr_patient_barcode <- 
                    make_bcr_code(dataSubt$patient_barcode)
            dataSubt$gender <- NULL
            dataSubt$age_at_initial_pathologic_diagnosis <- NULL
            clin_tumor_all <- merge(clin_tumor, dataSubt, 
                                    by = "bcr_patient_barcode",
                                    all = TRUE)
            } else {
            clin_tumor_all <- clin_tumor
            }
            clin_tumor_all <- 
            clin_tumor_all[!duplicated(clin_tumor_all$bcr_patient_barcode),]
        return(clin_tumor_all)
}


combine_info_for_cancer <- function(tumor_clin_file_path,
                                    tumor_pd1_dir = NULL) {
                                    # ind_scores_dir = NULL, this should 
                                    # pathways = NULL,
                                    # cluster_file = NULL) {

        # Combine clinical and tumor data
        clin_tumor <- combine_clins(tumor_clin_file_path)
        data <- list("clin" = clin_tumor)

        # Load PD1 data if pd1_dir is provided
        if (!is.null(tumor_pd1_dir)) {
            pd1_links <- load_pd1_generic(tumor_pd1_dir, type = "pd1_links")
            pd1_net <- load_pd1_generic(tumor_pd1_dir, type = "pd1_net")
            rownames(pd1_net) <- pd1_links
            data$pd1_net <- pd1_net
            data$pd1_scores <- load_pd1_generic(tumor_pd1_dir,
                                type = "pd1_scores")
            data$cd274_scores <- load_pd1_generic(tumor_pd1_dir, 
                                type = "cd274_scores")
            data$cd274_expression <- load_pd1_generic(tumor_pd1_dir, 
                                type = "cd274_expression")
            colnames(data$cd274_expression) <- gsub("\\."  , "-", 
                                colnames(data$cd274_expression))
        } else {
            message("Warning: 'tumor_pd1_dir' missing. PD1 data not included.")
        }

        # # Load pathway scores if both ind_scores_dir and pathways are provided
        # if (!is.null(ind_scores_dir) && !is.null(pathways)) {
        #     data$pathways_scores <- preprocess_pathways_scores(
        #         tumor, ind_scores_dir, pathways)
        # } else {
        #     if (is.null(ind_scores_dir)) 
        #         message("Warning: 'ind_scores_dir' missing.
        #                                     Pathways not included.")
        #     if (is.null(pathways))
        #         message("Warning: 'pathways' missing. 
        #                                     Pathways not included.")
        # }
        # if (!is.null(cluster_file)) {
        #         clusters <- read_in_cluster_info(cluster_file)
        #         clusters <- clean_cluster_info(clusters, tumor)
        #         clusters <- rearrange_cluster_info(clusters)
        #         data$clusters <- clusters
        #     } else {
        #         if (is.null(cluster_file))
        #         message("Warning: 'cluster_file' missing. 
        #                                     Clusters not included.")
        #     }

        return(data)
}

#' Prepare Data for a Glmnet Cox Model
#'
#' This function transforms feature data and combines it with clinical data to 
#' produce a single data frame suitable for fitting a Cox model.
#'
#' @param data_all A list containing:
#'   - `data`: A matrix where columns are individuals and rows are features.
#'   - `clin_data`: A data frame with clinical data for the individuals.
#'
#' @return A combined data of features with clinical data.
#'
#' @export

create_cox_data <- function(data_all) {
        data <- data_all$data
        clin_data <- data_all$clin_data
        features_t <- t(data_all$data)
        data_merged <- cbind(features_t, data_all$clin_data)
        return(data_merged)
}

#' Combine the clinical data with other data 
#' @param data cols are individuals, rows are features
#' @param clin_data should contain bcr_patient_barcode
#' @return 
#' @export 
#' 
#' 
create_data_combined <- function(data, clin_data) {
    olap <- which(colnames(data) %in% clin_data$bcr_patient_barcode)
    olap <- colnames(data)[olap]
    data <- data[, colnames(data) %in% olap]
    clin_data <- clin_data[clin_data$bcr_patient_barcode %in% olap, ]
    clin_data <- clin_data[match(colnames(data), 
                clin_data$bcr_patient_barcode), ]
    data_combined <- list("data" = data, "clin_data" = clin_data)
    return(data_combined)
}
#' Run multivariate Cox Model with LASSO Regularization
#'
#' Fits a Cox model to survival data using LASSO penalty with cross-validation.
#' Returns features with nonzero coefficients, indicating selected predictors.
#'
#' @param data_cox A data frame with `ToF_death` (time) and `event` (status).
#' @param number_folds Number of folds for cross-validation. Defaults to 5.
#' @param ntimes Number of iterations for stabilizing feature selection. 
#'               Defaults to 10.
#' @return A data table of selected features and their nonzero coefficients.

run_glmnet_cox_pd1 <- function(data_cox,
                                number_folds = 5,
                                ntimes = 10, alpha = 0.5){
        if (!all(c("ToF_death", "event") %in% colnames(data_cox))) {
            stop("`data_cox` must contain `ToF_death` and `event` columns.")
            }
        X <- as.matrix(data_cox[, grep("CD274", colnames(data_cox)), with = FALSE])
        if (ncol(X) == 0) {
            stop("No columns matching 'CD274' found in `data_cox`")
            }  
        res_all <- foreach::foreach(i = 1:ntimes, .combine = rbind) %dopar% {
            set.seed(i)
            foldid <- sample(rep(seq(number_folds), length.out = nrow(X)))
            fit_cox <- 
                glmnet::cv.glmnet(X, Surv(data_cox$ToF_death, data_cox$event),
                                family = "cox", 
                                standardize = T, 
                                nfolds =  number_folds, 
                                foldid = foldid,
                                alpha = alpha)
            lambda_optimal <- fit_cox$lambda.min
            coef_cox <- coef(fit_cox, s = lambda_optimal)
            coefficients <- coef_cox[coef_cox[,1] != 0, ]
            selected_edges <- rownames(coef_cox)[coef_cox[, 1] != 0]
            data_res <- data.table::data.table("edges" = selected_edges,
                                   "coefficients" = coefficients)
        }
        return(res_all)
}

#' Clean Data for Cox Regression Analysis
#'
#' Removes rows with non-positive `ToF_death` or missing `event`/`ToF_death`.
#'
#' @param data_cox A data frame with survival data. `ToF_death` is time to event
#' and `event` is the event indicator (1 = occurred).
#'
#' @return A data frame without rows with non-positive `ToF_death` or missing
#' values in `event`/`ToF_death`.
#' @export


clean_data_cox <- function(data_cox) {
    data_cox <- data_cox %>%
                dplyr::filter(ToF_death > 0) %>%
                dplyr::filter(!is.na(event) & !is.na(ToF_death))
}



run_glmnet_ntimes_pd1 <- function(
                            tumor_clin_file_path,
                            tumor_pd1_dir,
                            # ind_scores_dir = NULL,
                            # pathways = NULL, 
                            number_folds = 5,
                            ntimes = 10,
                            ncores = 10,
                            alpha = 1) {
    data <- combine_info_for_cancer(tumor_clin_file_path, 
                                tumor_pd1_dir)
                                # ind_scores_dir,
                                # pathways)
    data_combined <- create_data_combined(data$pd1_net, data$clin)
    data_cox <- create_cox_data(data_combined)
    doParallel::registerDoParallel(ncores)
    types <- c("OS", "PFI")
    selected_genes <- NULL
    selected_genes_all <- NULL
    for (i in 1:length(types)){
        type <- types[i]
        event_col <- paste0(type, collapse = "") 
        time_col <- paste0(type, ".time", collapse = "")
        data_cox$event <- as.numeric(data_cox[[event_col]])
        data_cox$ToF_death <- as.numeric(data_cox[[time_col]])
        data_cox <- clean_data_cox(data_cox)
        if (nrow(data_cox) == 0) next  # should be warning here
        if (nrow(data_cox[data_cox$event==1,])<5) next
        selected_genes <- 
                run_glmnet_cox_pd1(data_cox, number_folds, ntimes, alpha)
        selected_genes$type <- type
        selected_genes_all <- rbind(selected_genes, selected_genes_all)
    }
    return(selected_genes_all)
}