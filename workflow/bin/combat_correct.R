#####################
## Load R packages ##
#####################
required_libraries <- c("sva", "data.table", "optparse")
for (lib in required_libraries) {
    suppressPackageStartupMessages(
        library(lib, character.only = TRUE, quietly = TRUE)
    )
}

####################
## Read arguments ##
####################
option_list <- list(
    optparse::make_option(
        c("-e", "--expression_file"),
        type = "character",
        default = NULL,
        help = "Path to the normalized expression file.",
        metavar = "character"
    ),
    optparse::make_option(
        c("-g", "--group_file"),
        type = "character",
        default = NULL,
        help = "Path to the sample group file.",
        metavar = "character"
    ),
    optparse::make_option(
        c("-b", "--batch_file"),
        type = "character",
        default = NULL,
        help = "Path to the batch information file.",
        metavar = "character"
    ),
    optparse::make_option(
        c("-c", "--clin_file"),
        type = "character",
        default = NULL,
        help = "Path to the clinical file.",
        metavar = "character"),
    optparse::make_option(
        c("-o", "--output_file"),
        type = "character",
        default = NULL,
        help = "Path to the output file.",
        metavar = "character"
    )
)
opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

EXPRESSION_FILE <- opt$expression_file
GROUP_FILE <- opt$group_file
BATCH_FILE <- opt$batch_file
CLIN_FILE <- opt$clin_file
OUTPUT_FILE <- opt$output_file

source("workflow/bin/analyze_batch_fn.R")

cat("Loading expression file (this may take a while)...\n")
log2exp <- log2_transform_snail(EXPRESSION_FILE)

cat("Loading group file...\n")
groups <- fread(GROUP_FILE, head = FALSE)

cat("Loading batch file...\n")
batch_info <- fread(BATCH_FILE, head = TRUE)
batch_info <- batch_info[batch_info$Tissues %in% c("cancer")]
batch_info$Tissues <- NULL

cat("Loading clinical file...\n")
clin <- load_clin_rdata(CLIN_FILE)

batch_info$platform <- 
    clin$gdc_platform[match(
        batch_info$Samples, 
        clin$gdc_cases.samples.portions.analytes.aliquots.submitter_id
    )]
batch_info[is.na(batch_info)] <- c("not_available")

# Define cancer-specific configurations
cancer_configs <- list(
    "TCGA-COAD" = list(
        plates_to_remove = 2066, 
        apply_combat = TRUE, 
        verbose = FALSE
    ),
    "TCGA-DLBC" = list(
        plates_to_remove = 2404, 
        apply_combat = TRUE, 
        verbose = FALSE
    ),
    "TCGA-LUAD" = list(
        plates_to_remove = NULL, 
        apply_combat = TRUE, 
        verbose = FALSE
    ),
    "TCGA-PRAD" = list(
        plates_to_remove = 2302, 
        apply_combat = FALSE, 
        verbose = FALSE
    ),
    "TCGA-READ" = list(
        plates_to_remove = "A32Y", 
        apply_combat = TRUE, 
        verbose = FALSE
    ),
    "TCGA-UCEC" = list(
        plates_to_remove = NULL,
        apply_combat = TRUE, 
        verbose = FALSE
    )
)

# Get cancer types present in your data
available_cancers <- unique(groups$V2)
cat("Cancer types found in data:", 
    paste(available_cancers, collapse = ", "), "\n")

# Filter configs to only include cancers present in your data
cancers_to_correct <- intersect(names(cancer_configs), available_cancers)
cat("Cancer types that will be batch corrected:", 
    paste(cancers_to_correct, collapse = ", "), "\n")

# Handle case where no cancers need correction
if (length(cancers_to_correct) == 0) {
    cat("No cancer types found that require batch correction.\n")
    cat("Saving original data without correction...\n")
    
    # Filter to include only cancer samples
    log2exp_filtered <- log2exp[, colnames(log2exp) %in% batch_info$Samples]
    log2exp_new <- log2exp_filtered
    
} else {
    # Apply batch correction to cancers that are present and need correction
    corrected_data <- list()
    
    for (cancer in cancers_to_correct) {
        cat("Processing cancer type:", cancer, "\n")
        config <- cancer_configs[[cancer]]
        
        # Add error handling for each cancer type
        tryCatch({
            corrected_data[[cancer]] <- correct_cancer_batch(
                cancer_type = cancer,
                log2exp = log2exp,
                groups = groups,
                batch_info = batch_info,
                plates_to_remove = config$plates_to_remove,
                apply_combat = config$apply_combat,
                verbose = TRUE
            )
            cat("Successfully processed:", cancer, "\n")
        }, error = function(e) {
            cat("Warning: Failed to process", cancer, 
                "- Error:", e$message, "\n")
            cat("Skipping batch correction for", cancer, "\n")
        })
    }
    
    # Handle case where corrected_data might be empty
    if (length(corrected_data) > 0) {
        corrected_matrices <- as.data.frame(
            do.call("cbind", corrected_data)
        )
        cat("Corrected data dimensions:", dim(corrected_matrices), "\n")
        
        # Remove original data for cancer types successfully corrected
        cancers_successfully_corrected <- names(corrected_data)
        data_to_replace <- groups[
            groups$V2 %in% cancers_successfully_corrected, 
        ]
        log2exp_filtered <- log2exp[
            , !colnames(log2exp) %in% data_to_replace$V1
        ]
        
        # Include only cancer samples (remove healthy tissues)
        log2exp_filtered <- log2exp_filtered[
            , colnames(log2exp_filtered) %in% batch_info$Samples
        ]
        
        # Combine filtered data with corrected data
        cat("Combining corrected data with remaining samples...\n")
        log2exp_new <- cbind(log2exp_filtered, corrected_matrices)
        
    } else {
        cat("No cancer types were successfully corrected.\n")
        cat("Saving original data without correction...\n")
        
        # Filter to include only cancer samples
        log2exp_filtered <- log2exp[
            , colnames(log2exp) %in% batch_info$Samples
        ]
        log2exp_new <- log2exp_filtered
    }
}

cat("Final data dimensions:", dim(log2exp_new), "\n")
cat("Saving to:", OUTPUT_FILE, "\n")
save(log2exp_new, file = OUTPUT_FILE)
cat("Done!\n")