#' Log2 Transform SNAIL Normalized Expression Data
#'
#' Loads an RData file with SNAIL-normalized gene expression data,
#' rounds values to the nearest integer, and applies log2(x + 1).
#'
#' @param exp_file Path to the RData file with "exp" object. The first column is
#'   assumed non-numeric (e.g., gene IDs), others are numeric expression values.
#'
#' @return Matrix of log2-transformed expression values 
#'

log2_transform_snail <- function(exp_file) {
    exp <- fread(exp_file, header = F)
    colnms <- exp[2, ]
    exp <- exp[-c(1:2), ]
    colnames(exp) <- as.character(colnms)
    colnames(exp)[1] <- c("gene_id")
    exp_mat <- apply(exp[, -1], 2, function(x) as.numeric(x))
    exp_mat <- round(exp_mat)
    log2exp <- apply(exp_mat, 2, function(x) log2(x + 1))
    return(log2exp)
}

#' Load and clean batch information from a file
#'
#' Reads a batch info file and replaces NA in 'Year', 'Plates', 'TSS', and
#' 'Center' columns with "uknown".
#'
#' @param batch_file Path to the batch information file.
#'
#' @return Data frame with batch info and NAs replaced by "uknown".
#'
#' @import data.table
#' @export
load_batch <- function(batch_file) {
    batch_info <- fread(batch_file)
    batch_info$Year[is.na(batch_info$Year)] <- c("uknown")
    batch_info$Plates[is.na(batch_info$Plates)] <- c("uknown")
    batch_info$TSS[is.na(batch_info$TSS)] <- c("uknown")
    batch_info$Center[is.na(batch_info$Center)] <- c("uknown")
    return(batch_info)
}

#' Load Clinical Data from RData File
#'
#' Loads a clinical data object named \code{clin} from a specified RData file.
#'
#' @param clin_file Path to the RData file containing the \code{clin} object.
#'
#' @return The clinical data object \code{clin} loaded from the file.
#'
#' @export
load_clin_rdata <- function(clin_file) {
    clin_fpath <- file.path(clin_file)
    load(clin_fpath, clin <- new.env())
    clin <- clin[["clin"]]
    return(clin)
}





#' Batch Process Cancer Data
#'
#' Processes gene expression data for a cancer type, performing PCA and batch
#' effect analysis.
#'
#' @param tumor_type Character. Cancer type to process (e.g., "BRCA").
#' @param log2exp Matrix/data.frame. Log2-transformed expression (genes x samples).
#' @param groups Data.frame. Sample group info, columns V1 (IDs), V2 (project).
#' @param batch_info Data.frame. Batch info, must have 'Samples' column.
#' @param clin Data.frame. Clinical data (unused).
#' @param nthreads Integer. Threads for permutation tests (default: 1).
#' @param dsc_permutations Integer. DSC permutations for batch effect (default: 1000).
#' @param min_batch_size Integer. Minimum batch size (default: 2).
#' @param max_gene_count Integer. Max genes to include (default: 70000).
#' @param data_version Character. Data version (default: "1.0").
#' @param test_version Character. Test version (default: "1.0").

batch_process_cancer <- function(tumor_type,
                                log2exp,
                                groups, 
                                batch_info, 
                                nthreads = 1,
                                dsc_permutations = 1000,
                                min_batch_size = 2,
                                max_gene_count = 70000,
                                data_version = "1.0",
                                test_version = "1.0",
                                random_seed = 314) {
    cat(sprintf("Processing cancer type: %s\n", tumor_type))
    cat(sprintf("Current working directory: %s\n", getwd()))
    result <- tryCatch({
        proj <- paste("TCGA", tumor_type, sep = "-")
        samples <- groups$V1[groups$V2 == proj]
        # Check if samples exist for this cancer type
        if (length(samples) == 0) {
            cat(sprintf(
                "No samples found for %s, skipping...\n",
                tumor_type
            ))
            return(list(
                cancer = tumor_type,
                status = "no_samples",
                error = NULL
            ))
        }
        # Filter expression data for the current cancer type
        exp_proj <- log2exp[, colnames(log2exp) %in% samples, drop = FALSE]
        batches <- batch_info[batch_info$Samples %in% colnames(exp_proj), ]
        exp_proj <- 
            exp_proj[, match(batches$Samples, colnames(exp_proj)), drop = FALSE]
        # Check if the samples in expression data match the batch info
        if (!all(colnames(exp_proj) == batches$Samples)) {
            cat(sprintf(
                "Sample mismatch for %s, skipping...\n",
                tumor_type
            ))
            return(list(
                cancer = tumor_type,
                status = "sample_mismatch",
                error = NULL
            ))
        }
        # Check if there are enough samples for PCA
        # Minimum sample check
        if (ncol(exp_proj) < 10) {  # Minimum sample check
            cat(sprintf(
                "Too few samples for %s (%d), skipping...\n",
                tumor_type,
                ncol(exp_proj)
            ))
            return(list(
                cancer = tumor_type,
                status = "too_few_samples",
                error = NULL
            ))
        }

        myData <- new("BEA_DATA", as.matrix(exp_proj), batches, data.frame())
        PCA_Regular_Structures(theData = myData,
                             theTitle = paste('PCA', tumor_type, sep = '_'),
                             theOutputDir = ".",
                             theBatchTypeAndValuePairsToRemove = NULL,
                             theBatchTypeAndValuePairsToKeep = NULL,
                             theListOfComponentsToPlot = c(1, 2),
                             theDoDSCFlag = TRUE,
                             theDoDscPermsFileFlag = TRUE,
                             theDSCPermutations = dsc_permutations,
                             theDSCThreads = nthreads,
                             theMinBatchSize = min_batch_size,
                             theDataVersion = data_version,
                             theTestVersion = test_version,
                             theSeed = random_seed,
                             theMaxGeneCount = max_gene_count)
        # Log success  
        cat(sprintf("Successfully processed %s\n", tumor_type))
        return(list(
            cancer = tumor_type,
            status = "success",
            error = NULL
        ))
    }, error = function(e) {
        cat(sprintf("Error processing %s: %s\n", tumor_type, e$message))
        return(list(
            cancer = tumor_type,
            status = "error",
            error = e$message
        ))
    })
    return(result)
}
