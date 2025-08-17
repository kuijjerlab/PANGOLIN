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


#' Combine DSC results from multiple batches for a cancer type.
#'
#' Reads PCA annotation files from batch directories in output_dir/cancer.
#' Extracts DSC values and p-values for each batch, combines into a data.table.
#'
#' @param output_dir Path to output directory with batch results.
#' @param cancer Name of cancer type subdirectory to process.
#'
#' @return data.table with DSC values and p-values for each batch.
#'
#' @import data.table
#' @importFrom data.table fread
combine_mbatch_results <- function(output_dir, cancer) {
    files <- list.files(file.path(output_dir, cancer))
    files <- files[!grepl("ALL__CompListDSC.RData|Purity_singscore", files)]
    dsc_data_all <- NULL
    for (batch in files) {
        annotation_file <- file.path(
            output_dir, cancer, batch, 
            "ManyToMany", "1.0", "1.0", "PCAAnnotations.tsv"
        )
        if (!file.exists(annotation_file)) {
            warning(paste("File not found:", annotation_file))
            next
        }
        full_data <- fread(annotation_file)
        dsc_data <- data.table(
            "Annotation" = c("DSC all", "DSC (1,2)"),
            "DSC" = c(
                full_data[Annotation == "Disp. Sep. Crit. (DSC)", Value],
                full_data[Annotation == "Disp. Sep. Crit. (DSC) (1,2)", Value]
            ),
            "p_value" = c(
                full_data[Annotation == "DSC pvalue", Value],
                full_data[Annotation == "DSC pvalue(1,2)", Value]
            )
        )
        dsc_data$batch <- batch
        dsc_data$cancer <- cancer
        dsc_data_all <- rbind(dsc_data, dsc_data_all)
    }
    return(dsc_data_all)
}


#' Plot DSC values for batches/cancers with faceting and pagination.
#'
#' Creates a bar plot of DSC values for batches, faceted by cancer type and
#' paginated. Horizontal lines show DSC thresholds. Optionally saves to PDF.
#'
#' @param data Data frame with 'batch', 'DSC', and 'cancer' columns.
#' @param output_file Optional PDF file path for saving the plot.
#' @param colors Optional vector of fill colors for batches.
#' @param dsc_thresholds Numeric vector of DSC threshold lines.
#' @param plot_width Plot width (in inches).
#' @param plot_height Plot height (in inches).
#' @param facet_ncol Facet grid columns.
#' @param facet_nrow Facet grid rows.
#' @param page_num Facet page number.
#' @param text_size Axis text size.
#' @param line_size Threshold line thickness.
#' @param show_legend Logical, show legend.
#' @param angle_x_text X-axis text angle.
#'
#' @return ggplot object for DSC bar plot.
#'
#' @examples
#' # plot_mbatch_dsc(data = my_data, output_file = "dsc_plot.pdf")
#'
#' @import ggplot2
#' @importFrom ggforce facet_wrap_paginate
plot_mbatch_dsc <- function(
                data,
                output_file = NULL,
                colors = NULL,
                dsc_thresholds = c(0.6, 1.0),
                plot_width = 8,
                plot_height = 8,
                facet_ncol = 6,
                facet_nrow = 6,
                page_num = 1,
                text_size = 10,
                line_size = 0.5,
                show_legend = FALSE,
                angle_x_text = 90
            ) {
    if (is.null(colors)) {
        colors <- c(
            "#c85735", "#7c5fcd", "#83a23e", "#cc53ad", "#4eaa76",
            "#d33f63", "#5792ce", "#c28c42", "#a679bf", "#c06a7e"
        )
    }
    data$DSC <- as.numeric(data$DSC)
    p <- ggplot(
        data,
        aes(x = batch, y = DSC, fill = batch)
    ) +
        geom_bar(stat = "identity") +
        theme_bw() +
        scale_fill_manual(values = colors) +
        theme(
            axis.text.x = element_text(
                angle = angle_x_text,
                vjust = 0.5,
                hjust = 1,
                size = text_size
            ),
            legend.direction = "horizontal",
            legend.position = if (show_legend) "bottom" else "none"
        )
    for (threshold in dsc_thresholds) {
        p <- p + geom_hline(
            yintercept = threshold,
            linetype = "dashed",
            color = "black",
            size = line_size
        )
    }
    p <- p + facet_wrap_paginate(
        ~ cancer,
        ncol = facet_ncol,
        nrow = facet_nrow,
        page = page_num
    )
    if (!is.null(output_file)) {
        pdf(output_file, width = plot_width, height = plot_height)
        print(p)
        dev.off()
        cat("Plot saved to:", output_file, "\n")
    }
    return(p)
}

#' Perform batch correction using ComBat
#'
#' Applies ComBat from `sva` to correct batch effects in expression data.
#'
#' @param data Numeric matrix/data.frame, genes x samples.
#' @param batch_info Data.frame with batch info (e.g., year, sample IDs).
#' @param clin Data.frame with clinical info (platform, sample IDs).
#' @param batch_parameter Batch variable ("platform" or other).
#'
#' @return Batch-corrected expression matrix.
#' @importFrom sva ComBat
#' @export
combat_correct <- function(data, batch_info, clin, batch_parameter) {
    samples <- colnames(data)
    if (batch_parameter %in% c("platform")) {
        batch <- clin$gdc_platform[
            match(samples,
                  clin$gdc_cases.samples.portions.analytes.aliquots.submitter_id)
        ]
        cat("Samples without platform info:",
            length(batch[is.na(batch)]), "\n")
        print(samples[is.na(batch)])
        data <- data[, !is.na(batch)]
        batch <- batch[!is.na(batch)]
        batch <- as.factor(batch)
        combat_exp <- sva::ComBat(dat = data, batch = batch, mod = NULL,
                                  par.prior = TRUE, prior.plots = FALSE)
    } else {
        batch <- batch_info$Year[
            match(samples, batch_info$Samples)
        ]
        cat("Samples without year info:",
            length(batch[is.na(batch)]), "\n")
        print(samples[is.na(batch)])
        data <- data[, !is.na(batch)]
        batch <- batch[!is.na(batch)]
        batch <- as.factor(batch)
        combat_exp <- sva::ComBat(dat = data, batch = batch, mod = NULL,
                                  par.prior = TRUE, prior.plots = FALSE)
    }
    return(combat_exp)
}
