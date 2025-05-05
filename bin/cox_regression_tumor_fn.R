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

#' Extract PD-1 Pathway Individual Scores
#'
#' Loads an RData file containing individual scores and extracts
#' the PD-1 signaling pathway scores.
#'
#' @param individual_scores_file Path to the RData file with scores.
#' @return A data frame with PD-1 pathway scores, columns "PC1" & "PC2".
#' @export
extract_pd1_pathway_individual_scores <- function(individual_scores_file) {
        # Check if file exists
        if (!file.exists(individual_scores_file)) {
        stop("Error: File does not exist.")
        }
        data_env <- new.env()
        # Load the RData file
        tryCatch({
        load(individual_scores_file, envir = data_env)
        }, error = function(e) {
        stop("Error loading file: ", e$message)
        })
        # Validate 'ind_scores' presence
        if (!"ind_scores" %in% ls(data_env)) {
        stop("Error: 'ind_scores' not found.")
        }
        ind_scores <- data_env[["ind_scores"]]
        # Validate PD-1 signaling pathway presence
        if (!"REACTOME_PD_1_SIGNALING" %in% names(ind_scores)) {
        stop("Error: 'REACTOME_PD_1_SIGNALING' missing.")
        }
        # Extract and validate scores
        ind_scores_df <- as.data.frame(ind_scores$REACTOME_PD_1_SIGNALING)
        if (ncol(ind_scores_df) < 2) {
        stop("Error: Insufficient columns in data.")
        }
        rownames(ind_scores_df) <- make_bcr_code(rownames(ind_scores_df))
        colnames(ind_scores_df) <- c("PC1", "PC2")
        return(ind_scores_df)
        }

#' Load Curated Tumor Clinical Data
#'
#' This function reads a curated tumor clinical data file and returns it as a
#' data.table. 
#'
#' @param tumor_clin_file_path Character string specifying the path to the 
#' tumor clinical data file.
#'
#' @return A `data.table` containing the clinical data from the specified file.
#' @import data.table
#' @export

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

#' Load PD1 Generic Data
#'
#' This function loads and processes PD1 pathway data from a specified 
#' directory based on the provided type.
#'
#' @param tumor_pd1_dir Character string specifying the directory containing 
#' the PD1 data files.
#' @param type Character string indicating the type of PD1 data to load.
#'
#' @return A processed dataset containing PD1 data.
#' @export


load_pd1_generic <- function(tumor_pd1_dir, type) {
        validate_inputs(tumor_pd1_dir, type)
        file_pattern <- determine_pattern(type)
        tumor_file <- find_tumor_file(tumor_pd1_dir, file_pattern)
        data <- load_process_pd1_data(tumor_pd1_dir, tumor_file, type)
        return(data)
}

#' Validate Inputs for PD1 Data Loading
#'
#' Ensures that the provided directory exists and that the type argument is 
#' one of the allowed values.
#'
#' @param tumor_pd1_dir Character string specifying the directory to check.
#' @param type Character string indicating the type of PD1 data.
#'
#' @return No return value. Throws an error if validation fails.
#'
#' @export

validate_inputs <- function(tumor_pd1_dir, type) {
        if (!dir.exists(tumor_pd1_dir)) {
                stop("`pd1_dir` doesn't exist.")
        }
        if (!type %in% c("pd1_links", "pd1_net", 
                         "pd1_scores", "pdl1_expression")) {
                stop("Invalid `type`. Choose valid type.")
        }
}

#' Determine the Corresponding Pattern Name
#'
#' Maps a given `type` to its corresponding pattern from a predefined list.
#'
#' @param type A string. One of `"pd1_links"`, `"pd1_net"`, `"pd1_scores"`.
#' @return A string with the corresponding pattern or `NULL` if `type` is invalid.
#'
#' @examples
#' determine_pattern("pd1_links")  # Returns "pd1_edges"
#' determine_pattern("pd1_net")    # Returns "pd1_net"
#' determine_pattern("pd1_scores") # Returns "pd1_individual_scores"
#' determine_pattern("pdl1_expression")    # Returns "pdl1_expression"
#' determine_pattern("invalid")    # Returns NULL
#'
#' @export
#' 
determine_pattern <- function(type) {
        patterns <- list(pd1_links = "pd1_edges", pd1_net = "pd1_net", 
                         pd1_scores = "pd1_individual_scores",
                         pdl1_expression = "pdl1_expression")
        return(patterns[[type]])
}

#' Find a File Matching a Given Pattern
#'
#' Searches for a file in `tumor_pd1_dir` that matches the given `pattern`.
#' If no files or multiple files are found, an error is raised.
#'
#' @param tumor_pd1_dir A string. Directory path to search for the file.
#' @param pattern A string. Regular expression pattern to match filenames.
#' @return A string. The name of the matched file.
#' @throws Error if no file or multiple files match the pattern.
#' @export

find_tumor_file <- function(tumor_pd1_dir, pattern) {
        pattern_files <- list.files(tumor_pd1_dir, pattern = pattern, 
                                full.names = FALSE)
        if (length(pattern_files) == 0) {
                stop("No files found in: ", tumor_pd1_dir)
        }

        if (length(pattern_files) == 0) {
                stop("No file found in: ", tumor_pd1_dir)
        }
        if (length(pattern_files) > 1) {
                stop("Several files found: ", tumor_pd1_dir)
        }
        return(pattern_files[1])
}

#' Load and Process PD1 Data
#'
#' Loads a specified PD1 data file and processes it based on the given `type`.
#'
#' @param tumor_pd1_dir A string. Directory containing the PD1 data file.
#' @param tumor_file A string. Filename of the PD1 data file.
#' @param type A string. One of `"pd1_links"`, `"pd1_net"`, or `"pd1_scores"`.
#' @return A processed data object based on the selected `type`.
#' @throws Error if the file cannot be loaded or `type` is invalid.
#' @export
#' 
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
        } else if (type == "pdl1_expression") {
                data <- data.table::fread(file_path)
                }
        return(data)
}

#' Load an Object from an RData File
#'
#' Loads a specified object from an `.RData` file into a new environment.
#'
#' @param file A string. Path to the `.RData` file.
#' @param object_name A string. Name of the object to retrieve from the file.
#' @return The requested object from the `.RData` file.
#' @throws Error if the file does not exist or the object is not found.
#' @export
#' 
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
#' This function combines clinical information for a tumor type
#' 
#' @param tumor_clin_file_path Path to the clinical data file
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

#' Combine Clinical and Tumor Data for Analysis
#'
#' Merges clinical and tumor data, including PD-1 data if 
#' `tumor_pd1_dir` is provided.
#'
#' @param tumor_clin_file_path A string. Path to the clinical tumor data file.
#' @param tumor_pd1_dir A string (optional). Directory containing PD-1 data.
#' @return A list containing combined clinical and tumor data, including 
#' PD-1 data.

combine_info_for_cancer <- function(tumor_clin_file_path,
                                    tumor_pd1_dir = NULL,
                                    # ind_scores_dir = NULL, this should 
                                    # pathways = NULL,
                                    cluster_file = NULL) {

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
            data$pdl1_expression <- load_pd1_generic(tumor_pd1_dir,
                                type = "pdl1_expression")
        
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
        if (!is.null(cluster_file)) {
                clusters <- read_in_cluster_info(cluster_file)
                clusters <- clean_cluster_info(clusters, tumor)
                clusters <- rearrange_cluster_info(clusters)
                data$clusters <- clusters
            } else {
                if (is.null(cluster_file))
                message("Warning: 'cluster_file' missing. 
                                            Clusters not included.")
            }

        return(data)
}

#' Read Cluster Information
#'
#' This function reads a cluster file from cola clustering and 
#' processes its data.
#'
#' @param cluster_file Path to the cluster file.
#' @return A data frame with cluster information.
#' @import data.table
#' @export
read_in_cluster_info <- function(cluster_file) {
        # Check if file exists
        if (!file.exists(cluster_file)) {
                stop("Cluster file does not exist.")
        }

        # Attempt to read file
        cluster_df <- tryCatch({
                data.table::fread(cluster_file)
        }, error = function(e) {
                stop("Error reading cluster file: ", e$message)
        })

        # Process cluster information
        cluster_df$cluster <- paste0("cluster_", cluster_df$class)
        cluster_df$bcr_patient_barcode <- tryCatch({
                make_bcr_code(cluster_df$ID)
        }, error = function(e) {
                stop("Error in generating BCR code: ", e$message)
        })

        return(cluster_df)
}

#' Clean Cluster Information
#'
#' Filters cluster data for a specific tumor type.
#'
#' @param cluster_df Data frame with cluster information.
#' @param tumor Tumor type to filter by.
#' @return A cleaned data frame containing data for the specified tumor.
#' @export
#' 
clean_cluster_info <- function(cluster_df, tumor) {
        # Check input types
        if (!is.character(tumor) || length(tumor) != 1) {
                stop("Input tumor should be a single character string.")
        }
        tumor <-  suppressWarnings(toupper(tumor))
        # Filter by tumor type
        cluster_df_clean <- cluster_df[cluster_df$cancer == tumor, ]

        return(cluster_df_clean)
}

#' Rearrange Cluster Information
#'
#' Reshapes cluster information for analysis by rearranging data based on 
#' patient barcode.
#'
#' @param cluster_df Data frame with cluster information.
#' @return A matrix with clusters arranged by patient barcode.
#' @import data.table
#' @export
#' 
rearrange_cluster_info <- function(cluster_df) {
        # Reshape data
        clusters <- tryCatch({
                dcast(cluster_df, k ~ bcr_patient_barcode,
                            value.var = "cluster")
        }, error = function(e) {
                stop("Error in reshaping data: ", e$message)
        })

        # Store cluster IDs and remove column
        cluster_ids <- clusters$k
        clusters$k <- NULL
        clusters <- as.data.frame(clusters)
        rownames(clusters) <- cluster_ids
        return(clusters)
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
#' @export List of features and clinical data combined
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


#' Run GLMNET Cox Regression Multiple Times for PD-1 Data
#'
#' Performs Cox regression using GLMNET on PD-1 data, running the model `ntimes`
#' with cross-validation.
#'
#' @param tumor_clin_file_path A string. Path to the clinical tumor data file.
#' @param tumor_pd1_dir A string. Directory containing PD-1 network data.
#' @param number_folds An integer. Number of cross-validation folds
#' (default: 5).
#' @param ntimes An integer. Number of times to repeat the model
#' (default: 10).
#' @param ncores An integer. Number of CPU cores for parallel computing 
#' (default: 10).
#' @param alpha A numeric. Alpha Parameter for GLMNET (default: 1).
#' @return A data frame containing selected genes and 
#' their associated survival type.
#' @throws Warning if data is empty or contains fewer than 5 events.
#'
#' @export
run_glmnet_ntimes_pd1 <- function(
                        tumor_clin_file_path,
                        tumor_pd1_dir,
                        number_folds = 5,
                        ntimes = 10,
                        ncores = 10,
                        alpha = 1) {
    data <- combine_info_for_cancer(tumor_clin_file_path, 
                                tumor_pd1_dir)
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


#' Run Univariate Cox Proportional Hazards Model
#'
#' Performs univariate Cox regression on a selected data type, using specified 
#' covariates for survival analysis.
#'
#' @param tumor_clin_file_path A string. Path to the clinical tumor data file.
#' @param tumor_pd1_dir A string. Directory containing PD-1 data.
#' @param covariates A character vector. Covariates to include in the model.
#' @param datatype A string. Data type to analyze, default: `"pd1_scores"`.
#' @param type_stat A string. Type of statistic to return, 
#' default: `"full_stats"`.
#' @return A named list of Cox models for each feature 
#' in the selected data type.
#' @export
#' 
run_univariate_coxph_model <- function(tumor_clin_file_path,
                                tumor_pd1_dir,
                                covariates,
                                # ind_scores_dir = NULL,
                                # pathways = NULL,
                                cluster_file = NULL,
                                datatype = c("pd1_scores",
#                               "pd1_net")
                               "clusters"),
                                type_stat = c("full_stats",
                                "overall_pval")) {
        data <- combine_info_for_cancer(tumor_clin_file_path, 
                                tumor_pd1_dir, cluster_file)
        data_combined <- create_data_combined(data[[datatype]], data$clin)
        data_cox <- create_cox_data(data_combined)
        rownames(data_cox) <- data_cox$bcr_patient_barcode
        if (!is.null(covariates)) {
                rows_to_include <- rownames(na.omit(data_cox[covariates]))
                data_cox <- data_cox[rows_to_include, ]
        } else {
                data_cox <- data_cox
        }
        if (!is.null(covariates)) {
            exclude_vars <-
                sapply(data_cox[, covariates, drop = FALSE],
                    function(x) length(unique(x)) == 1)
            covariates <- covariates[!exclude_vars]
        }
        features <- rownames(data[[datatype]])
        cox_models <- 
            purrr::map(features,
             ~create_univariate_cox_model(.x, covariates, data_cox, type_stat)) 
        names(cox_models) <- features
        return(cox_models)
}

#' Create Univariate Cox Proportional Hazards Model for a Given Feature
#'
#' @description Runs a Cox model for a feature with covariates and 
#' survival data.
#'
#' @param feature Character; feature for which the Cox model is created.
#' @param covariates Character; covariate(s) to include in the model.
#' @param data_cox Data frame; contains time-to-event and covariate information.
#'
#' @return List of Cox model summaries for OS (Overall Survival) 
#'          and PFI (Progression-Free Interval).


create_univariate_cox_model <- function(feature,
                                        covariates,
                                        data_cox,
                                        type_stat = c("full_stats",
                                                    "overall_pval")) {
        formulaString <- 
            paste("Surv(ToF_death, event) ~",
                    paste(covariates, collapse = " + "))
        feature <- sprintf("`%s`", feature)
        # Create the formula
        formula <- as.formula(paste(formulaString, "+", feature))
        coxph_model_list <- list()
        predicted_risk <- NULL
        types <- c("OS", "PFI")
        for (i in 1:length(types)){
            type <- types[i]
            event_col <- paste0(type, collapse = "")
            time_col <- paste0(type, ".time", collapse = "")
            data_cox$event <- as.numeric(data_cox[[event_col]])
            data_cox$ToF_death <- as.numeric(data_cox[[time_col]])
            data_cox <- clean_data_cox(data_cox)
            if (nrow(data_cox) == 0) next
            m <- "run without warning"
            coxph_model <- tryCatch(
                    {
                        survival::coxph(formula, data = data_cox)
                    },
                    warning = function(w) {
                        # Print the warning message for inspection
                        message("Warning occurred: ", conditionMessage(w))
                        m <<- "with warning"
                        survival::coxph(formula, data = data_cox)
                    }
                    )
            predicted_risk <- 
                        predict(coxph_model, newdata = data_cox, type = "risk")
            coxph_model_summary <- unpack_summary_coxph(coxph_model,
                                                type_stat = type_stat)
            coxph_model_summary$type <- type
            coxph_model_summary$m <- m
            predicted_risk <- as.data.frame(predicted_risk)
            predicted_risk$bcr_patient_barcode <- rownames(predicted_risk)
            predicted_risk$type <- type
            coxph_model_list[[i]] <- list("coxph_model" = coxph_model_summary,
                                    "predicted_risk" = predicted_risk)

        }
        return(coxph_model_list)
}

#' Unpack Summary of Cox Proportional Hazards Model
#'
#' This function extracts key statistics from a fitted Cox Proportional Hazards 
#' model (`coxph` model) and returns either the full statistics for each 
#' feature or only the overall p-value, depending on the specified output type.
#'
#' @param coxph_model A `coxph` model object created by the `coxph` function 
#'        from the `survival` package. This model should represent a Cox 
#'        Proportional Hazards regression fit.
#' @param type_stat A character string indicating the type of statistics to 
#'        return. Options are `"full_stats"` to return a table with hazard 
#'        ratios, p-values, and confidence intervals for each feature, or 
#'        `"overall_pval"` to return only the overall p-value for the model. 
#'        Default is `"full_stats"`.
#'
#' @return A `data.table` object. If `type_stat` is `"full_stats"`, the 
#'         returned table contains columns `feature`, `pvalue`, `hr` (hazard 
#'         ratio), `hr_lower` (lower confidence bound), and `hr_upper` (upper 
#'         confidence bound) for each feature. If `type_stat` is 
#'         `"overall_pval"`, the table contains a single column `pval`, which 
#'         is the overall p-value for the model.
#'
unpack_summary_coxph <- function(coxph_model, 
                                 type_stat = c("full_stats", "overall_pval")) {
        # Validate input arguments
        if (missing(coxph_model) || class(coxph_model) != "coxph") {
                stop("coxph_model must be a valid 'coxph' model object.")
        }
        type_stat <- match.arg(type_stat)

        # Extract summary details from coxph model
        coxph_model_summary <- summary(coxph_model)
        pvalue_res <- coxph_model_summary$coefficients[, "Pr(>|z|)"]
        hr_lower <- signif(coxph_model_summary$conf.int[, "lower .95"], 5)
        hr_upper <- signif(coxph_model_summary$conf.int[, "upper .95"], 5)
        hr <- exp(coef(coxph_model))
        overall_pvalue <- coxph_model_summary$logtest['pvalue']

        # Output based on the selected type_stat
        if (type_stat == "full_stats") {
                res <- data.table::data.table(
                        feature = names(hr),
                        pvalue = pvalue_res,
                        hr = hr,
                        hr_lower = hr_lower,
                        hr_upper = hr_upper
                )
        } else if (type_stat == "overall_pval") {
                res <- data.table::data.table(pval = overall_pvalue)
        }

        return(res)
}

#' Process PD-1 Univariate Cox Regression Results
#'
#' Extracts and combines Cox regression model results and predicted risk 
#' data from PD-1 analysis.
#'
#' @param pd1_res A list. Contains Cox regression results for PD-1 data.
#' @param pc_names A character vector. Principal component names to process.
#' @return A list with two data frames:
#'   \item{coxph_model_data}{Combined Cox model results.}
#'   \item{predicted_risk_data}{Combined predicted risk data.}
#' @export
#' 
process_pd1_univarite_cox_res <- function(pd1_res, pc_names) {
        for (pc in pc_names) {
        for (i in seq_along(pd1_res[[pc]])) {
        pd1_res[[pc]][[i]]$coxph_model$component <- pc
        pd1_res[[pc]][[i]]$predicted_risk$component <- pc
        }
        }
        # Combine coxph_model data
        coxph_model_data <- do.call(rbind, lapply(pc_names, function(pc) {
        do.call(rbind, lapply(pd1_res[[pc]], function(x) x$coxph_model))
        }))

        # Combine predicted_risk data
        predicted_risk_data <- do.call(rbind, lapply(pc_names, function(pc) {
        do.call(rbind, lapply(pd1_res[[pc]], function(x) x$predicted_risk))
        }))

        return(list(
        coxph_model_data = coxph_model_data,
        predicted_risk_data = predicted_risk_data
        ))
        }

#' Extract Gene Expression for a Specific Gene
#'
#' This function extracts the expression values of a specific gene
#'
#' @param exp_file Path to the gene expression file (CSV/TSV format). 
#'   The file should have genes in the first column and expression values 
#'   in subsequent columns.
#' @param samples_file Path to the samples file containing 
#'   sample IDs. 
#' @param gene_id The ID of the gene to extract (default: "CD274").
#' @return A data frame with the extracted gene expression values
#' @examples
#' @export
extract_gene_expression <- function(exp_file,
                        samples_file, 
                        gene_id = "CD274") {
                            if (!file.exists(exp_file)) {
        stop("The specified `exp_file` does not exist.")
        }
        if (!file.exists(samples_file)) {
        stop("The specified `samples_file` does not exist.")
        }
        if (!is.character(gene_id) || length(gene_id) != 1) {
        stop("`gene_id` must be a single string.")
        }
        gene_id_mod <- paste0("^", gene_id, "$")
        exp <- fread(exp_file)
        genes <- exp$V1
        exp <- exp[,-1]
        samples <- fread(samples_file)
        colnames(exp) <- samples$sample_id
        idx_gene <- grep(gene_id_mod, genes)
        if (length(idx_gene) == 0) {
        stop(paste("Gene ID", gene_id, "not found in the expression file."))
        }
        gene_exp <- t(exp[idx_gene,])
        colnames(gene_exp) <- gene_id
        gene_exp <- as.data.frame(gene_exp)
        gene_exp$cancer <- 
                samples$cancer[match(rownames(gene_exp), samples$sample_id)]
        gene_exp$bcr_patient_barcode <- make_bcr_code(rownames(gene_exp))
        return(gene_exp)
        }

