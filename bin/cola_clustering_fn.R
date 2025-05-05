set.seed(1234)
#' Perform COLA clustering on cancer data
#'
#' This function performs consensus clustering on cancer data, using either
#' "indegree" or "expression" data based on the specified type.
#'
#' @param cancer A character string specifying the cancer type to analyze.
#' @param exp_file A character string with the file path to the expression data.
#'        Only used when `datatype` is set to "expression".
#' @param indegree_dir A character string with the directory for indegree data.
#'        Only used when `datatype` is set to "indegree".
#' @param datatype A character string, either "indegree" or "expression",
#'        specifying the data type for clustering. Default is "indegree".
#' @param n_cores An integer for the number of cores to use for clustering.
#'        Default is 10.
#'
#' @return A result object from `consensus_clustering` with clustering details.
#' @export
#' 
perform_cola_clustering <- function(cancer,
                            exp_file = NULL,
                            samples_file = NULL,
                            indegree_dir = NULL,
                            datatype = c("indegree", "expression"),
                            n_cores = 1,
                            top_value_method = "ATC",
                            partition_method = "kmeans",
                            max_k = 6,
                            p_sampling = 0.8,
                            partition_repeat = 1000,
                            scale_rows = TRUE,
                            output_dir) {
    datatype <- match.arg(datatype)
    if (datatype == "indegree") {
        if (is.null(indegree_dir)) {
            stop("indegree directory must be provided.")
        }
        # List files and filter by cancer type
        ind_file <- list.files(indegree_dir, full.names = TRUE) %>%
            .[stringr::str_detect(., regex(cancer, ignore_case = TRUE))]

        if (length(ind_file) == 0) {
            stop("No indegree file found")
        }
        # Load indegree data
        ind_data <- tryCatch({
            load_indegree(ind_file)
        }, error = function(e) {
            stop("Error loading indegree data: ", e$message)
        })
        data <- ind_data[, -1]
    } else if (datatype == "expression") {
        if (is.null(exp_file)) {
            stop("exp_file must be provided for expression data.")
        }
        if (is.null(samples_file)) {
            stop("samples_file must be provided for expression data.")
        }
        # Load expression data
        exp_data <- tryCatch({
            load_exp(cancer, exp_file, samples_file)
        }, error = function(e) {
            stop("Error loading expression data: ", e$message)
        })
        data <- exp_data[, -1]
    }
        res <- perform_consensus_clustering(data, 
                                            n_cores = n_cores,
                                            top_value_method = top_value_method,
                                            partition_method = partition_method,
                                            max_k = max_k,
                                            p_sampling = p_sampling,
                                            partition_repeat = partition_repeat,
                                            scale_rows = TRUE)
    cat("saving results and plots", "\n")
    if (!dir.exists(output_dir)) {
        dir.create(output_dir, recursive = TRUE)
    }
        # Save results and generate plots
    save_results(res, cancer, datatype, output_dir)
}


#' Perform Consensus Clustering on Data
#'
#' This function performs consensus clustering using the COLA package on the
#' provided data matrix.
#'
#' @param data A data frame or matrix with samples in columns and features in rows.
#' @param n_cores Integer specifying the number of cores for parallel processing.
#'        Default is 10.
#' @param top_value_method Method for selecting top features, default is "ATC".
#' @param partition_method Clustering method, default is "hclust".
#' @param max_k Integer specifying the maximum number of clusters, default is 6.
#' @param p_sampling Proportion of samples used in each partitioning iteration.
#'        Default is 0.8.
#' @param partition_repeat Number of partitioning repetitions, default is 1000.
#' @return A consensus clustering result object.
#' @export
perform_consensus_clustering <- function(data,
                                         n_cores = 10,
                                         top_value_method = "ATC",
                                         partition_method = c("kmeans"),
                                         max_k = 6,
                                         p_sampling = 0.8,
                                         partition_repeat = 1000,
                                         scale_rows = TRUE) {
    # Input validation
    if (is.null(data) || ncol(data) < 2) {
        stop("Input data must be a matrix with at least two columns.")
    }
    mad <- as.matrix(data)
    res <- tryCatch({
        cola::consensus_partition(mad,
                                cores = n_cores,
                                top_value_method = top_value_method,
                                partition_method = partition_method,
                                max_k = max_k,
                                p_sampling = p_sampling,
                                partition_repeat = partition_repeat,
                                scale_rows = TRUE)
    }, error = function(e) {
        stop("Consensus clustering failed: ", e$message)
    })
    return(res)
}

#' Save Results and Generate Plots for Consensus Clustering
#'
#' This function saves results and generates relevant plots for a consensus 
#' clustering analysis. It creates multiple output files, including a summary 
#' of selected clusters, membership, statistics, and various plots, all named 
#' consistently based on the specified cancer type and data type.
#'
#' @param res A clustering result object from consensus clustering.
#' @param cancer A character string specifying the cancer type, used in file names.
#' @param datatype A character string, either "indegree" or "expression", 
#'        indicating the type of data used. Default is "indegree".
#' 
#' @return None. Saves various output files to disk.
#'
#' @export
save_results <- function(res, 
                cancer,
                datatype = c("indegree", "expression"),
                output_dir) {
        # Ensure datatype is valid
        datatype <- match.arg(datatype)
        # Determine the best and optional cluster values for k
        best_k <- suggest_best_k(res)
        selected_k <- best_k[1]
        optional_k <- if (!is.null(attr(best_k, "optional"))) 
                        attr(best_k, "optional") else NA

        # Helper function to create consistent file name prefix
        get_filename <- function(suffix) {
            file.path(output_dir, 
                paste0(cancer, "_", suffix, "_", datatype, ".RData"))
        }
        # 1. Save best k information
        res_k <- data.table(
            cancer = cancer,
            selected_k = selected_k,
            optional_k = optional_k
        )
        tryCatch({
            save(res_k, file = get_filename("best_k"))
        }, error = function(e) {
            warning("Failed to save best_k information: ", e$message)
        })
        # 2. Save the main results object
        tryCatch({
            save(res, file = get_filename("results"))
        }, error = function(e) {
            warning("Failed to save main results object: ", e$message)
        })
        # 3. Generate and save collected plots PDF
        tryCatch({
            pdf(file = file.path(output_dir, paste0(cancer, "_collected_plots_", datatype, ".pdf")),
                width = 10, height = 10)
            collect_plots(res)
            dev.off()
        }, error = function(e) {
            warning("Failed to generate collected plots PDF: ", e$message)
        })

        # 4. Generate and save t-SNE plot PDF
        tryCatch({
            pdf(file = file.path(output_dir, 
                paste0(cancer, "_tSNE_", datatype, ".pdf")),
                width = 8, height = 8)
            dimension_reduction(res, selected_k)
            dev.off()
        }, error = function(e) {
            warning("Failed to generate t-SNE plot PDF: ", e$message)
        })
        # 5. Save membership information
        tryCatch({
            membership <- get_membership(res, selected_k)
            save(membership, file = get_filename("membership"))
        }, error = function(e) {
            warning("Failed to save membership information: ", e$message)
        })
        # 6. Save statistics
        tryCatch({
            statistics <- get_stats(res)
            save(statistics, file = get_filename("statistics"))
        }, error = function(e) {
            warning("Failed to save statistics: ", e$message)
        })
        # 7. Generate and save partition selection plot PDF
        tryCatch({
            pdf(file = file.path(output_dir, 
                paste0(cancer, "_select_partition_", datatype, ".pdf")),
                width = 8, height = 8)
            select_partition_number(res)
            dev.off()
        }, error = function(e) {
            warning("Failed to generate partition selection: ", e$message)
        })
        # 8. Save classes information
        tryCatch({
            classes <- get_classes(res, k = selected_k)
            save(classes, file = get_filename("classes"))
        }, error = function(e) {
            warning("Failed to save classes information: ", e$message)
        })
    }

#' Load a specified object from an .RData file
#'
#' This function loads an .RData file into a new environment and retrieves
#' a specified object. 
#'
#' @param result_file A character string representing the path to the .RData file.
#' @param object_name A character string for the name of the object to load 
#'        from the file. Default is "res_k".
#' @return The specified object from the .RData file.


load_result <- function(result_file, object_name = "res_k") {
        # Check if file exists
        if (!file.exists(result_file)) {
                stop("The specified file does not exist: ", result_file)
        } 
        # Load the .RData file into a new environment
        env <- new.env()
        load(result_file, envir = env)
        # Check if the specified object exists in the environment
        if (!exists(object_name, envir = env)) {
                stop("The object '", object_name, "' was not found in the file.")
        }
        # Retrieve and return the object
        result <- env[[object_name]]
        return(result)
}

#' Extract relevant classes from result files for specified tumors
#'
#' This function loads and processes data from specified result files 
#' for multiple tumor types.
#'
#' @param cola_dir A character string specifying the directory containing 
#'        result files.
#' @param best_k_df_file A character string representing the path to a 
#'        file containing the best cluster information for each tumor.
#' @return A data frame containing extracted class assignments for 
#'         each tumor type.

extract_relevant_classes <- function(res_files, best_k_df_file) {
        data_long <- load_and_prepare_data(best_k_df_file)
        # Error if no result files are found
        if (length(res_files) == 0) {
                stop("Error: No result files found in cola_dir.")
        }
        # Process and extract classes for each tumor type
        extracted_classes <- process_all_tumors(data_long, res_files)
        return(extracted_classes)
}


#' Load and prepare clustering data
#'
#' This function loads a data file containing information on the best 
#' clusters for various tumor types.
#'
#' @param best_k_df_file A character string representing the path to the 
#'        data file containing cluster information for each tumor type.
#' @return A data frame where each row represents a single tumor and possible 
#'         number of clusters
#' 
load_and_prepare_data <- function(best_k_df_file) {
        # Check if file exists
        if (!file.exists(best_k_df_file)) {
                stop("Error: File does not exist. Check the file path.")
        }
        # Load data and prepare for processing
        data <- fread(best_k_df_file)
        data_long <- data %>%
                separate_rows(combined, sep = ";") %>%
                mutate(possible_clusters = as.integer(combined))
        return(data_long)
}

#' Process and extract classes for all tumors
#'
#' This function iterates over each row in the input data, processes 
#' each tumor using a list of result files, and extracts relevant class 
#' assignments.
#'
#' @param data_long A data frame where each row represents a tumor and 
#'        its associated information, including possible clusters.
#' @param res_files A character vector of file paths to result files for 
#'        different tumors.
#' @return A data frame containing extracted class assignments for all 
#'         tumors.
#'
process_all_tumors <- function(data_long, res_files) {
        extracted_classes <- NULL
        # Loop over each row in the data
        for (i in 1:nrow(data_long)) {
                data_long_i <- data_long[i, ]
                tumor_classes <- process_single_tumor(data_long_i, res_files)
                # Append the extracted classes
                extracted_classes <- rbind(tumor_classes, extracted_classes)
        }
        return(extracted_classes)
}


#' Process and extract classes for a single tumor
#'
#' This function processes a single tumor type by loading its associated 
#' result file, retrieving class assignments based on the best cluster 
#' (`best_k`), and returning the extracted classes.
#'
#' @param data_long_i A data frame row containing information for a 
#'        specific tumor, including the cancer type and possible clusters.
#' @param res_files A character vector of file paths to result files 
#'        for various tumors.
#' @return A data frame with extracted class assignments for the tumor,
#'         or NULL if the result file or class data is not found.
process_single_tumor <- function(data_long_i, res_files) {
        tumor <- data_long_i$cancer
        best_k <- data_long_i$possible_clusters

        # Find the relevant result file
        res_file <- res_files[grep(tumor, res_files)]
        if (length(res_file) == 0) {
                warning(paste("Warning: No result file found for tumor type", 
                              tumor))
                return(NULL)
        }
        # Load the result object
        res <- load_result(res_file, object_name = "res")
        # Retrieve class assignments based on best_k
        classes <- get_classes(res, best_k)
        if (is.null(classes)) {
                warning(paste("Warning: Failed to retrieve classes for tumor", 
                              "type", tumor, "and k =", best_k))
                return(NULL)
        }
        classes$k <- paste("k", best_k, sep = "_")
        classes$cancer <- tumor
        classes$ID <- rownames(classes)
        return(classes)
}


#' Perform t-SNE analysis for a specific tumor type
#'
#' This function locates the result file for a specified tumor, loads 
#' the necessary data, and performs a t-SNE analysis using PCA.
#'
#' @param tumor A character string representing the tumor type to analyze.
#' @param res_files A character vector of file paths to result files for 
#'        different tumors.
#' @param perplexity An integer specifying the perplexity parameter for t-SNE.
#'        Default is 10.
#' @param n_pcs An integer specifying the number of principal components to 
#'        use for dimensionality reduction before t-SNE. Default is 10.
#' @param scale A logical value indicating whether to scale the data before
#'        running t-SNE. Default is TRUE.
#' @return A data frame or matrix containing the t-SNE results for the 
#'         specified tumor.
#' @throws A warning if no or multiple result files are found for the 
#'         tumor, or if data is missing in the loaded result object.
#' @examples
#' # Example usage:
#' # perform_tsne_cancer("lung_cancer", res_files, perplexity = 30, n_pcs = 50)
#'
perform_tsne_cancer <- function(tumor, res_files, perplexity = 10, 
                                n_pcs = 10, scale = TRUE) {
        # Find the result file for the specified tumor
        res_file <- res_files[grep(tumor, res_files)]
        # Check for errors if no or multiple files are found
        if (length(res_file) == 0) {
                warning("No result file found for tumor type: ", tumor)
                return(NULL)
        } else if (length(res_file) > 1) {
                warning("Multiple result files found for tumor type: ", tumor, 
                        ". Using the first match.")
                res_file <- res_file[1]  # Use the first matching file
        }
        # Load the result object from the file
        res <- tryCatch({
                load_result(res_file, object_name = "res")
        }, error = function(e) {
                warning("Failed to load result for tumor type: ", tumor, 
                        " - ", e$message)
                return(NULL)
        })
        # Extract data and perform t-SNE analysis
        data <- res@.env$data
        tsne_res <- tryCatch({
                runTSNE_withPCA(cancer = tumor, data = data, 
                                perplexity = perplexity, 
                                n_pcs = n_pcs,
                                scale = scale)
        }, error = function(e) {
                warning("t-SNE analysis failed for tumor type: ", tumor, 
                        " - ", e$message)
                return(NULL)
        })
        return(tsne_res)
}

#' Perform t-SNE analysis for multiple tumor types
#'
#' This function iterates over a list of tumor types, locates their result 
#' files, and performs a t-SNE analysis for each tumor using data from the 
#' specified directory.
#'
#' @param tumors A character vector with the names of tumor types to analyze.
#' @param cola_dir A character string specifying the directory that contains 
#'        result files for different tumors.
#' @return A data frame combining the t-SNE results for all specified tumors.
#'         Returns NULL if no results are generated.
#' @throws An error if the 'tumors' list or 'res_files' list is empty.
#'         Warnings if t-SNE fails for any specific tumor.
#' @examples
#' # Example usage:
#' # perform_tsne_all_cancers(c("lung", "breast"), "path/to/cola_dir")
#'
perform_tsne_all_cancers <- function(tumors, res_files) {
        # Retrieve result files in the specified directory
        # res_files <- list.files(
        #         cola_dir, pattern = "results", full.names = TRUE
        # )
        # Check if tumors list and result files list are non-empty
        if (length(tumors) == 0) {
                stop("Error: The 'tumors' list is empty.")
        }
        if (length(res_files) == 0) {
                stop("Error: The 'res_files' list is empty.")
        }

        # Initialize container for all t-SNE results
        tsne_res_all <- NULL
        # Loop through each tumor and perform t-SNE
        for (i in seq_along(tumors)) {
                tumor <- tumors[i]
                print(tumor)

                # Try to perform t-SNE and handle errors
                tsne_res <- tryCatch({
                        perform_tsne_cancer(tumor, res_files)
                }, error = function(e) {
                        warning(paste("Warning: Failed to perform t-SNE for", 
                                      tumor, "-", e$message))
                        return(NULL)
                })

                # Accumulate results if t-SNE was successful
                if (!is.null(tsne_res)) {
                        tsne_res_all <- rbind(tsne_res, tsne_res_all)
                }
        }
        return(tsne_res_all)
}


#' This function performs PCA on input data to reduce dimensionality, 
#' then applies t-SNE on the principal components for visualization.
#'
#' @param cancer A character string indicating the cancer type.
#' @param data A matrix or data frame where rows are genes/features and 
#'        columns are samples.
#' @param perplexity A numeric value specifying the t-SNE perplexity 
#'        parameter.
#' @param n_pcs An integer specifying the number of principal components 
#'        to retain for t-SNE.
#' @param scale Logical, whether to scale data before PCA. Default is TRUE.
#' @return A data frame with t-SNE dimensions, sample IDs, and cancer type.
#' @throws An error if input data is missing, not a matrix, or has 
#'         insufficient dimensions.
#' @examples
#' # Example usage:
#' # runTSNE_withPCA("lung_cancer", data, perplexity = 30, n_pcs = 10)
#'
runTSNE_withPCA <- function(cancer, data, perplexity, n_pcs, scale = TRUE) {
        # Check for valid input data
        if (is.null(data) || !is.matrix(data) && !is.data.frame(data)) {
                stop("Error: 'data' must be a matrix or data frame.")
        }
        if (ncol(data) < n_pcs) {
                stop("Error: 'data' has fewer columns than 'n_pcs'.")
        }
        # Transpose data for PCA and t-SNE, handle potential errors
        data_t <- tryCatch({
                t(data)
        }, error = function(e) {
                stop("Error transposing data: ", e$message)
        })
        # Perform PCA with error handling
        pca_res <- tryCatch({
                irlba::prcomp_irlba(data_t, n_pcs, center = TRUE, scale. = scale)
        }, error = function(e) {
                stop("Error performing PCA: ", e$message)
        })
        # Extract principal components
        PCA <- pca_res$x
        # Perform t-SNE with error handling
        tsne_res <- tryCatch({
                Rtsne::Rtsne(PCA[, 1:n_pcs], perplexity = perplexity)
        }, error = function(e) {
                stop("Error performing t-SNE: ", e$message)
        })
        # Format t-SNE results into a data frame
        tsneData <- data.frame(tsne_res$Y)
        colnames(tsneData) <- c("Dim1", "Dim2")
        tsneData$id <- colnames(data)
        tsneData$cancer <- cancer
        return(tsneData)
}

#' Load indegree file
#' @param indegree_file path to indegree file
#' @return indegree
#' @export 

load_indegree <- function(indegree_file) {
    load(indegree_file, indegree <- new.env())
    indegree <- indegree[['ind']]
    return(indegree)
    }

#' Load expression
#' @param cancer we are interested in
#' @param exp_file expression file
#' @param samples_file samples file
#' @return indegree
#' @export
load_exp <- function(cancer, exp_file, samples_file) {
    exp <- fread(exp_file)
    genes <- exp$V1
    exp <- exp[,-1]
    samples <- fread(samples_file)
    colnames(exp) <- samples$sample_id
    samples_cancer <- 
                samples$sample_id[grep(cancer, samples$cancer, ignore.case = T)]
    exp_cancer <- exp[, colnames(exp) %in% samples_cancer, with=FALSE]
    exp_cancer <- cbind("gene" = genes, exp_cancer)
    exp_cancer <- exp_cancer[rowSums(exp_cancer[,-1]) >0, ]
    return(exp_cancer)
    }

#' Combine Best and Optional Cluster Results
#'
#' This function combines the selected and optional cluster values (`selected_k`
#' and `optional_k`) for each cancer type into a single string, separated by 
#' semicolons. It ensures that the values are unique and sorted.
#'
#' @param data A `data.table` with columns `cancer`, `selected_k`, and 
#'        `optional_k`.
#' @return A `data.table` with one row per cancer type and a `combined` column 
#'         containing the combined cluster values as a semicolon-separated 
#'         string.
#' @export

combine_k_results <- function(data){
    combined <- data[, {
        all_vals <- c(na.omit(selected_k), na.omit(optional_k))
        .(combined = paste(sort(unique(all_vals)), collapse = ";"))
            }, by = cancer]
}


#' Plot t-SNE Results for Consensus Clusters
#'
#' This function generates a t-SNE plot for consensus clusters using ggplot2.
#'
#' @param res_data A data frame containing t-SNE results and cluster information. 
#'   It must include the columns `Dim1`, `Dim2`, `class`, `cancer`, and `k`.
#' @param output_file A string specifying the path to save the PDF file.
#' @param plot_width Numeric, the width of the output PDF.
#' @param plot_height Numeric, the height of the output PDF.
#' @return A ggplot object representing the t-SNE plot.
#' @examples
#' # Example usage:
#' plot_tsne_clusters(res_ind, "output.pdf", 10, 14)
#' @export

plot_tsne_clusters_cola <- function(res_data) {
        POINT_SIZE <- 0.4
        POINT_ALPHA <- 0.6
        GRID_LINEWIDTH <- 0.5
        GRID_COLOR <- "grey80"
        STRIP_TEXT_SIZE <- 12
        AXIS_TITLE_SIZE <- 14
        AXIS_TEXT_SIZE <- 7
        LEGEND_TEXT_SIZE <- 10
    # Create the plot
    p <- ggplot(res_data, aes(x = Dim1, y = Dim2, color = class)) +
        geom_point(size = POINT_SIZE, alpha = POINT_ALPHA) +
        scale_color_viridis(discrete = TRUE, option = "D",
                            guide = guide_legend(nrow = 1)) +
        facet_wrap(~cancer + k, scales = "free") +
        theme_minimal() +
        theme(
        panel.grid.major = 
            element_line(linewidth = GRID_LINEWIDTH, color = GRID_COLOR),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        strip.text = element_text(size = STRIP_TEXT_SIZE, face = "bold"),
        axis.title = element_text(size = AXIS_TITLE_SIZE),
        axis.text = element_text(size = AXIS_TEXT_SIZE),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = LEGEND_TEXT_SIZE),
        legend.box = "horizontal"
        ) +
        labs(
        title = "",
        x = "t-SNE Dimension 1",
        y = "t-SNE Dimension 2"
        )
    # Return the ggplot object
    return(p)
}


