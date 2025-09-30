#' Load Indegree Data from an RData File
#'
#' This function loads an RData file containing an "indegree" object
#' and extracts the stored variable.
#'
#' @param indegree_file A string specifying the path to the RData file.
#' @return The loaded "indegree" object.

load_indegree <- function(indegree_file) {
        if (!file.exists(indegree_file)) {
            stop("File does not exist: ", indegree_file)
        }
        indegree_env <- new.env()
        load(indegree_file, envir = indegree_env)
        if (!exists("ind", envir = indegree_env)) {
            stop("'ind' object not found in the loaded file.")
        }
        return(indegree_env$ind)
}



#' Combine Indegree Data from Multiple Files
#'
#' This function loads and combines indegree data from multiple files
#'  in a specified directory.
#'
#' @param cancer_dir A path to the directory containing indegree files.
#' @return A data frame with combined indegree data.
#' @importFrom dplyr bind_cols
#' 
#' @export

combine_indegree <- function(filenames) {
        # Check if any files were found
        if (length(filenames) == 0) {
            stop("No files matching 'indegree_norm' 
            were found in the specified directory.")
        }
        # Load and combine data
        data_list <- lapply(filenames, load_indegree)

        data <- dplyr::bind_cols(data_list)

        # Identify and remove duplicate "tar" columns, keeping the first one
        tar_cols <- grep("tar", colnames(data))
        if (length(tar_cols) > 1) {
        data <- data[, -tar_cols[-1]]
        }
        # Rename the first "tar" column if necessary
        colnames(data)[1] <- "tar"
        return(data)
}

#' Perform t-SNE on PCA-reduced data
#'
#' This function first applies PCA to reduce the dimensionality
#'  of the input data, and then performs t-Distributed Stochastic 
#' Neighbor Embedding (t-SNE) using the `Rtsne` package.
#'
#' @param data A numeric matrix or data frame with features
#'  as rows and samples as columns.
#' @param perplexity Numeric; perplexity parameter for t-SNE, 
#' typically between 5 and 50. Default: 20
#' @param n_pcs Integer; number of principal components to retain for t-SNE.
#' Default: 20
#'
#' @return A data frame containing t-SNE coordinates (`Dim1`, `Dim2`)
#'  and sample IDs (`id`).
#' @export
runTSNE_withPCA <- function(data, perplexity = 20, n_pcs = 20) {
        if (!is.matrix(data) && !is.data.frame(data)) {
            stop("Input 'data' must be a matrix or data frame.")
            }
        if (!is.numeric(perplexity) || perplexity <= 0) {
            stop("'perplexity' must be a positive numeric value.")
            }
        if (!is.numeric(n_pcs) || n_pcs <= 0) {
            stop("'n_pcs' must be a positive integer.")
            }
        # Transpose data to have samples as rows
        data_t <- t(as.matrix(data))
        df <- data_t
        samples <- colnames(data)
        # Perform PCA using irlba
        PCA <- irlba::prcomp_irlba(df, n_pcs, center = T, scale. = T)$x
        # Run t-SNE on PCA-reduced data
        tsne_result <- Rtsne::Rtsne(PCA[, 1:n_pcs], 
                    perplexity = perplexity)$Y
        tsneData <- data.frame("Dim1" = tsne_result[, 1],
                    "Dim2" = tsne_result[, 2],
                    "id" = samples)
        return(tsneData)
}
