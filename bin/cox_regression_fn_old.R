# old file
#########################################################################
## Functions performing survival analysis on PD1 network edges, PCs 
## from PD1 subnetwork, and PC from PDL1 subnetwork (CD274)
##########################################################################

# libraries <- c("data.table", "dplyr", "doParallel",
#             "survminer", "survival", "magrittr", "gridExtra",
#             "purrr")
# lapply(libraries, library, character.only = TRUE)


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

#' Preprocess Immune Data
#'
#' This function reads in immune data from a file (CIBERSORTx) and preprocesses it.
#'
#' @param immune_file A string representing the path to the immune data file.
#' @return A `data.table` containing the preprocessed immune data.
#' @import data.table
#' @export
preprocess_immune <- function(immune_file) {
        # Check if the file exists
        if (!file.exists(immune_file)) {
            stop("The provided immune file doesn't exist.")
        }
        immune_data <- tryCatch({
            data.table::fread(immune_file)
        }, error = function(e) {
            stop("Failed to read the file: ", e$message)
        })

        # Check if the required column `Mixture` exists
        if (!"Mixture" %in% colnames(immune_data)) {
            stop("The input file must contain a column named `Mixture`.")
        }
        # Generate the `sample_id_r` column using `make_bcr_code`
        immune_data$sample_id_r <- tryCatch({
            make_bcr_code(immune_data$Mixture)
        }, error = function(e) {
            stop("Error generating `sample_id_r` column: ", e$message)
        })
        immune_data[, c("P-value", "Correlation", "RMSE")] <- NULL
        return(immune_data)
    }

#' Load curated clinical data
#'
#' This function loads curated clinical data from an Excel file.
#' The data file is downloaded from the article.
#' available at https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6066282/.
#'
#' @param clin_file_path A string representing the path 
#'                       to the clinical data Excel file.
#' @return A data frame containing the clinical data.
#' @import openxlsx
#' @export
load_clin_curated <- function(clin_file_path) {
        # Validate input: Check if the file exists
        if (!file.exists(clin_file_path)) {
            stop("The specified file path does not exist: ", clin_file_path)
        }
        clin <- tryCatch({
            openxlsx::read.xlsx(clin_file_path, sheet = "TCGA-CDR")
        }, error = function(e) {
            stop("Error reading the Excel file: ", e$message)
        })
        clin <- clin[, -1]
        return(clin)
}


#' Select clinical data for a specific tumor type
#'
#' This function filters clinical data to include only info corresponding to a
#' specified tumor type. 
#'
#' @param clin A data frame containing clinical data and must include
#'             a column named `type` that specifies the tumor type for each row.
#' @param tumor A string representing the tumor type to filter by (e.g., "acc").
#' @return A data frame containing clinical data only for the specified tumor.
#' @export
#' 
select_tumor_clin_curated <- function(clin, tumor) {
        # Validate input: Check if `tumor` is a non-empty string
        if (!is.character(tumor) || length(tumor) != 1 || tumor == "") {
            stop("The `tumor` argument must be a non-empty string.")
        }
        # Filter clinical data by the specified tumor type
        clin_tumor <- clin[grepl(tumor, clin$type, ignore.case = TRUE), ]
        # Warn if no matches are found
        if (nrow(clin_tumor) == 0) {
            warning("No matches found for tumor type: ", tumor)
        }
        # Return the filtered clinical data
        return(clin_tumor)
    }

#' Combine clinical, PD1 network, PD1 individual mapping scores,  
#' CD274 individual mapping scores, and pathways scores for a specific tumor type.
#' 
#' This function loads and combines clinical data, PD1 network information,
#' PD1 scores, and CD274 scores for a specified tumor type.
#' 
#' @param tumor The tumor type to gather information for
#' @param clin Clinical data
#' @param pd1_dir The directory containing PD1 files
#' @param ind_scores_dir directory containing individual mappings PORCUPINE
#'                       results for all pathways
#' @param ind_scores_dir directory containing individual mappings PORCUPINE
#'                       results for all pathways
#' @param pathways pathways to include for the analysis
#' @param cluster_file file with cluster results
#'  
#' @return A list containing combined information:
#'         clinical data, PD1 network, PD1 scores, and CD274 scores
#' @export
#' 
combine_info_for_cancer <- function(tumor, clin,
                                    pd1_dir = NULL,
                                    ind_scores_dir = NULL,
                                    pathways = NULL,
                                    cluster_file = NULL) {
        # Check required parameters
        if (is.null(tumor)) stop("Error: 'tumor' data must be provided.")
        if (is.null(clin)) stop("Error: 'clin' data must be provided.")

        # Combine clinical and tumor data
        clin_tumor <- combine_clins(tumor, clin)
        data <- list("clin" = clin_tumor)

        # Load PD1 data if pd1_dir is provided
        if (!is.null(pd1_dir)) {
            pd1_links <- load_pd1_generic(tumor, pd1_dir, type = "pd1_links")
            pd1_net <- load_pd1_generic(tumor, pd1_dir, type = "pd1_net")
            rownames(pd1_net) <- pd1_links
            data$pd1_net <- pd1_net
            data$pd1_scores <- load_pd1_generic(tumor, pd1_dir,
                                type = "pd1_scores")
            data$cd274_scores <- load_pd1_generic(tumor, pd1_dir, 
                                type = "cd274_scores")
            data$cd274_expression <- load_pd1_generic(tumor, pd1_dir, 
                                type = "cd274_expression")
            colnames(data$cd274_expression) <- gsub("\\."  , "-", 
                                colnames(data$cd274_expression))
        } else {
            message("Warning: 'pd1_dir' missing. PD1 data not included.")
        }

        # Load pathway scores if both ind_scores_dir and pathways are provided
        if (!is.null(ind_scores_dir) && !is.null(pathways)) {
            data$pathways_scores <- preprocess_pathways_scores(
                tumor, ind_scores_dir, pathways)
        } else {
            if (is.null(ind_scores_dir)) 
                message("Warning: 'ind_scores_dir' missing.
                                            Pathways not included.")
            if (is.null(pathways))
                message("Warning: 'pathways' missing. 
                                            Pathways not included.")
        }
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
#' Load and Process the PD1 Data
#'
#' This function serves as a generic pipeline to load and process PD1 data 
#' based on the specified tumor, directory, and analysis type. It validates 
#' inputs, determines the appropriate subdirectory and file pattern, locates 
#' the tumor file, and processes the data accordingly.
#'
#' @param tumor A string specifying the tumor identifier. Must be a non-empty 
#'   string.
#' @param pd1_dir A string specifying the base directory containing PD1 data.
#' @param type A string specifying the type of PD1 analysis. Valid options are:
#'   `"pd1_links"`, `"pd1_net"`, `"pd1_scores"`, and `"cd274_scores"`.
#'
#' @return The function returns the processed PD1 data object:
#'   - For `"pd1_links"`, a vector of concatenated regulator-target pairs.
#'   - For `"pd1_net"`, a PD1 pathway subnetwork
#'   - For `"pd1_scores"` and `"cd274_scores"`, a matrix of mappings of
#'     individuals on the corresponding PCs based on pd1 pathway or 
#'     only the PDL1 gene
#'
#' @export

# load_pd1_generic <- function(tumor, pd1_dir, type) {
#         validate_inputs(tumor, pd1_dir, type)
#         pd1_subdir <- determine_subdir(pd1_dir, type)
#         file_pattern <- determine_pattern(type)
#         tumor_file <- find_tumor_file(pd1_subdir, tumor, file_pattern)
#         data <- load_process_pd1_data(pd1_subdir, tumor_file, type)
#         return(data)
# }

load_pd1_generic <- function(tumor, pd1_dir, type) {
        validate_inputs(tumor, pd1_dir, type)
        # pd1_subdir <- determine_subdir(pd1_dir, type)
        file_pattern <- determine_pattern(type)
        tumor_file <- find_tumor_file(pd1_dir, tumor, file_pattern)
        data <- load_process_pd1_data(pd1_dir, tumor_file, type)
        return(data)
}

#' Validate Input Parameters for PD1 Analysis
#'
#' This function validates the input parameters for PD1-related analysis.
#' It ensures the tumor identifier is a non-empty string, the directory
#' exists, and the `type` parameter is one of the valid options.
#'
#' @param tumor A non-empty string representing the tumor identifier.
#' @param pd1_dir A string for the path to the directory with PD1 data. Must
#'                exist.
#' @param type A string specifying the analysis type. Valid options are:
#'             `"pd1_links"`, `"pd1_net"`, `"pd1_scores"`, and "cd274_scores"`.
#'
#' @return Stops execution and raises an error if inputs are invalid.
#' @export
#' 
validate_inputs <- function(tumor, pd1_dir, type) {
        if (!is.character(tumor) || length(tumor) != 1) {
                stop("`tumor` must be a non-empty string.")
        }
        if (!dir.exists(pd1_dir)) {
                stop("`pd1_dir` doesn't exist.")
        }
        if (!type %in% c("pd1_links", "pd1_net", 
                         "pd1_scores", "cd274_scores",
                         "cd274_expression")) {
                stop("Invalid `type`. Choose valid type.")
        }
}

#' Determine Subdirectory for PD1 Analysis
#'
#' This function determines the subdirectory within the provided PD1 directory
#' based on the specified analysis type. It validates the existence of the
#' subdirectory and raises an error if it is not found.
#'
#' @param pd1_dir A string specifying the base directory for PD1 data.
#' @param type A string specifying the analysis type. Valid options are:
#'             `"pd1_links"`, `"pd1_net"`, `"pd1_scores"`, and
#'              `"cd274_scores"`.
#'
#' @return  A string representing the path to the subdirectory.
#'          Stops with an error if the subdirectory doesn't exist.
#'
#' @export
#' 
# determine_subdir <- function(pd1_dir, type) {
#         subdirs <- list(pd1_links = "pd1_edges",
#                         pd1_net = "pd1_edges",
#                         pd1_scores = "pd1_scores", 
#                         cd274_scores = "pd1_scores",
#                         cd274_expression = "pdl1_expression")
#         pd1_subdir <- file.path(pd1_dir, subdirs[[type]])
#         if (!dir.exists(pd1_subdir)) {
#                 stop("Subfolder not found: ", pd1_subdir)
#         }
#         return(pd1_subdir)
# }

#' Determine File Pattern for PD1 Analysis
#'
#' This function returns the file naming pattern associated with the specified
#' analysis type. The pattern corresponds to the expected file naming 
#' conventions used for each type of PD1 analysis.
#'
#' @param type A string specifying the analysis type. Valid options are:
#'   `"pd1_links"`, `"pd1_net"`, `"pd1_scores"`, and `"cd274_scores"`.
#'
#' @return A string representing the file naming pattern for the specified 
#'   analysis type
#'
#' @export
#' 
determine_pattern <- function(type) {
        patterns <- list(pd1_links = "pd1_edges", pd1_net = "pd1_net", 
                         pd1_scores = "pd1_individual_scores", 
                         cd274_scores = "CD274_individual_scores",
                         cd274_expression = "pdl1_expression")
        return(patterns[[type]])
}

#' Find Tumor File in PD1 Subdirectory
#'
#' This function searches for a file in a given PD1 subdirectory that matches
#' both the specified file naming pattern and tumor identifier. It ensures that
#' the file exists and raises an error if no matching file is found.
#'
#' @param pd1_subdir A string specifying the subdirectory path to search within.
#' @param tumor      A string specifying the tumor identifier to match in the
#'                   file name. Case-insensitive matching is applied.
#' @param pattern    A string specifying the file naming pattern to filter files
#'                   in the directory.
#'
#' @return A string representing the name of the first file that matches the
#'   tumor identifier and file naming pattern.
#' 
# find_tumor_file <- function(pd1_subdir, tumor, pattern) {
#         dir_files <- list.files(pd1_subdir, pattern = pattern, 
#                                 full.names = FALSE)
#         if (length(dir_files) == 0) {
#                 stop("No files found in: ", pd1_subdir)
#         }
#         tumor_file <- dir_files[grep(tumor, dir_files, ignore.case = TRUE)]
#         if (length(tumor_file) == 0) {
#                 stop("No file for tumor: ", tumor)
#         }
#         return(tumor_file[1])
# }


find_tumor_file <- function(pd1_dir, tumor, pattern) {
        dir_files <- list.files(pd1_dir, pattern = pattern, 
                                full.names = FALSE)
        if (length(dir_files) == 0) {
                stop("No files found in: ", pd1_dir)
        }
        tumor_file <- dir_files[grep(tumor, dir_files, ignore.case = TRUE)]
        if (length(tumor_file) == 0) {
                stop("No file for tumor: ", tumor)
        }
        return(tumor_file[1])
}


#' Load and Process PD1 Data
#'
#' This function loads and processes PD1 data from a file in a specified 
#' subdirectory. The processing steps depend on the analysis type, which 
#' determines the object to load and the specific transformations applied.
#'
#' @param pd1_subdir A string specifying the path to the PD1 subdirectory 
#'   containing the tumor file.
#' @param tumor_file A string specifying the name of the tumor file to load.
#' @param type A string specifying the type of PD1 analysis. Valid options are:
#'   `"pd1_links"`, `"pd1_net"`, `"pd1_scores"`, and `"cd274_scores"`.
#'
#' @return A processed data object, depending on the type of PD1 analysis:
#' 
#' @export
# load_process_pd1_data <- function(pd1_subdir, tumor_file, type) {
#         file_path <- file.path(pd1_subdir, tumor_file)
#         if (type == "pd1_links") {
#                 data <- load_pd1_object(file_path, object_name = "links")
#                 data <- paste(data$reg, data$tar, sep = "_")
#         } else if (type == "pd1_net") {
#                 data <- load_pd1_object(file_path, object_name = "pd1_net")
#                 colnames(data) <- make_bcr_code(colnames(data))
#         } else if (type == "pd1_scores") {
#                 scores <- load_pd1_object(file_path, object_name = "ind_scores")
#                 data <- t(scores)
#                 colnames(data) <- make_bcr_code(colnames(data))
#                 data <- data[1:2, ] # only take 2 PCs
#         } else if (type == "cd274_scores") {
#                 scores <- load_pd1_object(file_path, object_name = "ind_scores")
#                 data <- t(scores)
#                 colnames(data) <- make_bcr_code(colnames(data))
#         } else if (type == "cd274_expression") {
#                 data <- load_pd1_object(file_path, object_name = "pdl1_expression")
#                 data <- data.frame(data)
#         }
#         return(data)
# }

load_process_pd1_data <- function(pd1_dir, tumor_file, type) {
        file_path <- file.path(pd1_dir, tumor_file)
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

#' Load a Specific PD1 Object from an RData File
#'
#' This function loads a specified object from an RData file. It ensures the 
#' file exists and checks that the requested object is present in the file.
#'
#' @param file A string specifying the path to the RData file.
#' @param object_name A string specifying the name of the object to load from 
#'   the RData file.
#'
#' @return The function returns the requested object from the RData file. If
#'   the file does not exist or the object is not found, it raises an error.
#'
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
combine_clins <- function(tumor, clin) {
        if (!"bcr_patient_barcode" %in% colnames(clin)) {
            stop("The `clin` data frame must contain a column 
                named `bcr_patient_barcode`.")
        }
        if (!is.character(tumor) || length(tumor) != 1 || tumor == "") {
            stop("The `tumor` argument must be a non-empty string.")
        }
        tumor <-  suppressWarnings(toupper(tumor))
        clin_tumor <- tryCatch({select_tumor_clin_curated(clin, tumor)
            }, error = function(e) {
            stop("Error filtering clinical data for tumor type: ", e$message)
            })
        # those tumors are not included in TCGAquery_subtype
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
        X <- as.matrix(data_cox[, grep("CD274", colnames(data_cox))])
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

#' Run LASSO Regularized Cox Regression Multiple Times
#'
#' Prepares data for Cox regression, runs LASSO-regularized models repeatedly.
#'
#' @param tumor Tumor type (e.g., "acc").
#' @param clin A data frame with clinical data.
#' @param pd1_dir Directory with PD1 files.
#' @param ncores Number of cores for parallel processing. Defaults to 10.
#' @param number_folds Number of folds for cross-validation. Defaults to 5.
#' @param ntimes Times to run Cox model for stable feature selection. 
#'               Defaults to 10
#'
#' @return Data frame of selected genes with non-zero coefficients from all runs
#' of Cox models for Overall Survival and Progression-Free Interval.
#' @importFrom doParallel registerDoParallel
#' @export

run_glmnet_ntimes_pd1 <- function(tumor, 
                            clin,
                            pd1_dir,
                            ind_scores_dir = NULL,
                            pathways = NULL, 
                            number_folds = 5,
                            ntimes = 10,
                            ncores = 10,
                            alpha = 0.5) {
    data <- combine_info_for_cancer(tumor, 
                                clin, 
                                pd1_dir, 
                                ind_scores_dir,
                                pathways)
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
        selected_genes <- run_glmnet_cox_pd1(data_cox, number_folds, ntimes, alpha)
        selected_genes$type <- type
        selected_genes_all <- rbind(selected_genes, selected_genes_all)
    }
    return(selected_genes_all)
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


#' Run Univariate Cox Proportional Hazards Models for Each Feature
#'
#' This function prepares the dataset by combining tumor data, clinical data,
#'  and PD1-related data and then fits univariate Cox proportional hazards models 
#' for each feature (gene) while adjusting for specified covariates.
#'
#' @param tumor tumortype
#' @param clin A data frame containing clinical data.
#' @param pd1_dir A directory path or reference to PD1-related data.
#' @param covariates A character vector specifying the names of covariates 
#'                   to adjust for in the Cox models.
#' @param datatype A character vector of data types to use. 
#'                 The options are 'pd1_net', 'pd1_scores', 'cd274_scores', 
#'                 "pathways_scores".
#'                 The function will use the selected datatype to extract
#'                 feature data from the combined dataset.
#'
#' @return A named list of Cox proportional hazards models. Each model corresponds to a
#'         separate feature (gene) adjusted for specified covariates.
#'
#' @import purrr
#' @import dplyr 
#' @import survival 
#' @export
#'

run_univariate_coxph_model <- function(tumor,
                                clin,
                                pd1_dir,
                                covariates,
                                ind_scores_dir = NULL,
                                pathways = NULL,
                                cluster_file = NULL,
                                datatype = c("pd1_net",
                                "pd1_scores",
                                "cd274_scores",
                                "pathways_scores",
                                "clusters"),
                                type_stat = c("full_stats",
                                "overall_pval")) {
        data <- combine_info_for_cancer(tumor, clin, pd1_dir = pd1_dir, 
                ind_scores_dir = ind_scores_dir, pathways = pathways,
                cluster_file = cluster_file)
        data_combined <- create_data_combined(data[[datatype]], data$clin)
        data_cox <- create_cox_data(data_combined)
        rows_to_include <- rownames(na.omit(data_cox[covariates]))
        #data_cox$clin <- data_cox[rows_to_include, ]
        data_cox <- data_cox[rows_to_include, ]
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

#' Combine Data from a Nested List into a Single Data Table
#'
#' @description Combines data tables from a nested list into one, adding a 
#' `cluster` column.
#'
#' @param res Named list of lists, where each contains data tables for a cluster.
#'
#' @return A data table with all rows, including a `cluster` column indicating 
#' their origin.
#'
#' @importFrom data.table rbindlist
#'
#' @export
#' 
combine_data_from_list <- function(res) {
        flattened_list <- lapply(names(res), function(cluster) {
                data_list <- res[[cluster]]
                lapply(data_list, function(dt) {
                        dt[, cluster := cluster]
                        return(dt)
                })
        })
        flattened_list <- unlist(flattened_list, recursive = FALSE)
        final_df <- data.table::rbindlist(flattened_list,
                                    use.names = TRUE, fill = TRUE)
        return(final_df)
}


#' Analyze Cox model coefficients with LASSO regularization results.
#'
#' Aggregates LASSO results by TF, computes median coefficients, and counts TFs.
#'
#' @param res A data frame containing the model results.
#' @param threshold Minimum occurrences of a TF to include in the summary.
#'        Default: 1.
#'
#' @return A data table with three columns: 'TF' (factor name), 'coeff' (median
#'         coefficient), and 'ntimes'. Filter by 'threshold'.
#'
#' @import plyr
#' @import data.table
#' @export
#' 

analyze_res_coefficients <- function(res, threshold = 1) {
        res$TF <- sapply(strsplit(res$edges, "_"), function(x) x[1])
        coeff <- plyr::ddply(res, .(TF), function(x) median(x$coefficients))
        ntimes <- plyr::ddply(res, .(TF), function(x) nrow(x))
        summary_data <- data.table::data.table("TF" = coeff$TF, 
                                "coeff" = coeff$V1,
                                "ntimes" = ntimes$V1) 
        summary_data <- summary_data[summary_data$ntimes >= threshold, ]
        return(summary_data)

}

#' Read pathway scores for a specific tumor type.
#'
#' Reads pathway scores for a tumor type from the specified directory.
#'
#' @param tumor Tumor type (e.g., "acc").
#' @param ind_scores_dir Directory with individual scores files.
#' @return Data frame of pathway scores for the tumor type.
#' @export
#' 
read_pathways_scores <- function(tumor, ind_scores_dir) {
        # List files in the directory matching "individual_scores"
        ind_scores_files <- list.files(
            file.path(ind_scores_dir), pattern = "individual_scores"
        )
        # Filter files by tumor type, ignoring case
        ind_scores_file <- ind_scores_files[
            grep(tumor, ind_scores_files, ignore.case = TRUE)
        ]
        # Handle errors: No matching file or multiple matches
        if (length(ind_scores_file) == 0) {
            stop("No matching file found for tumor type: ", tumor)
        }
        load(file.path(ind_scores_dir, ind_scores_file), data <- new.env())
        if (!"ind_scores" %in% ls(data)) {
            stop("'ind_scores' not found in the loaded file.")
        }
        data <- data[["ind_scores"]]
        return(data)
    }

#' Filter Pathway Scores
#'
#' Filters pathway scores data to include only the specified pathways.
#'
#' @param data A data frame with pathway scores.
#' @param pathways Character vector of pathways to include.
#' @return Data frame with scores for specified pathways.
#' @export
#' 
filter_pathway_scores <- function(data, pathways) {
        data_clean <- data[names(data) %in% pathways]
        if (nrow(data_clean) == 0) {
            warning("No matching pathways found in the data.")
        }
        return(data_clean)
}


#' Read and Preprocess Pathway Scores
#'
#' Read and preprocess pathway scores for a given tumor type.
#'
#' @param tumor The tumor type.
#' @param ind_scores_dir Directory with individual scores files.
#' @param pathways A vector of pathways to filter in the data.
#' @export
#' 
preprocess_pathways_scores <- function(tumor,
                                        ind_scores_dir,
                                        pathways) {
    if (!is.character(tumor) || length(tumor) != 1) {
        stop("`tumor` must be a single string.")
        }
    if (!dir.exists(ind_scores_dir)) {
        stop("`ind_scores_dir` does not exist.")
        }
    if (!is.character(pathways) || length(pathways) < 1) {
        stop("`pathways` must be a non-empty character vector.")
        }
    pathways_scores <- read_pathways_scores(tumor, ind_scores_dir)
    pathways_scores <- 
        filter_pathway_scores(pathways_scores, pathways)
    pathways_scores <- lapply(names(pathways_scores), function(name) {
        pathway_scores <- pathways_scores[[name]]
        colnames(pathway_scores) <- 
            paste(name, colnames(pathway_scores), sep = "_")
        rownames(pathway_scores) <- make_bcr_code(rownames(pathway_scores))
        pathway_scores
    })
    pathways_scores <- t(do.call(cbind, pathways_scores))
    return(pathways_scores)
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
        gene_exp$cancer <- samples$cancer[match(rownames(gene_exp), samples$sample_id)]
        gene_exp$bcr_patient_barcode <- make_bcr_code(rownames(gene_exp))
        return(gene_exp)
        }

extract_gene_expression_tumor <- function(exp_file, 
                                           samples_file, 
                                           gene_id = "CD274", 
                                           tumor = tumor) {
    
        # Extract gene expression data
        gene_exp <- extract_gene_expression(exp_file, samples_file, gene_id)
        gene_exp$cancer <- gsub("TCGA-", "", gene_exp$cancer)
        # Ensure `cancer` column exists in `gene_exp`

        # Filter for the specific tumor type
        filtered_data <- gene_exp[gene_exp$cancer == tumor, ]
        
        # Pivot data to wide format
        wide_data <- filtered_data %>%
                pivot_wider(names_from = bcr_patient_barcode, values_from = gene_id) %>%
                mutate(gene = gene_id) %>%
                select(gene, everything())
        wide_data <- as.data.frame(wide_data[,-c(1:2)])
        rownames(wide_data) <- gene_id
        return(wide_data)
}

#' Create a Scatter Plot of Principal Component vs. CD274 Expression
#'
#' Generates a scatter plot showing the relationship between a principal 
#' component (`pc_component`) and CD274 (PD-L1) expression. The points 
#' are colored by `risk_score`, and a correlation statistic is displayed.
#'
#' @param data A data frame containing the relevant variables.
#' @param pc_component A string specifying the column name for the principal 
#'   component.
#' @param cancer A string representing the cancer type, used as the plot title.
#'
#' @return A ggplot2 object representing the scatter plot.
#' @export

create_pc_cd274_plot <- function(data, pc_component, cancer) {
    # Validate required columns
    required_columns <- c(pc_component, "CD274_exp", "risk_score")
    missing_cols <- setdiff(required_columns, names(data))
    if (length(missing_cols) > 0) {
        stop(paste("Error: Missing columns:", 
                paste(missing_cols, collapse = ", ")))
    }

    # Compute median for color scaling
    midpoint_value <- median(na.omit(data$risk_score))

    # Create scatter plot
    ggplot(data = data, aes(x = .data[[pc_component]], y = CD274_exp)) +
        geom_point(aes(color = risk_score), size = 2) +
        scale_color_gradient2(low = "blue", mid = "white", high = "red",
                              midpoint = midpoint_value) +  
        theme_minimal() +
        sm_statCorr(corr_method = "spearman", color = "black") +
        ylab("PDL1 (CD274) expression") +
        ggtitle(cancer)
}


# create_pc_cd274_plot <- function(data, pc_component, cancer) {
#         # Compute the correct median for color scaling
#         midpoint_value <- median(na.omit(data$risk_score))
#         ggplot(data = data, aes(x = !!sym(pc_component), y = CD274_exp)) +
#         geom_point(aes(color = risk_score), size = 2) +
#         scale_color_gradient2(low = "blue", mid = "white", high = "red",
#                         midpoint = median(midpoint_value)) +
#         theme_minimal() +
#         sm_statCorr(corr_method = "spearman", color = "black") +
#         ylab("PDL1 (CD274) expression") +
#         ggtitle(cancer)
# }
