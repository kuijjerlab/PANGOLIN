#####################
## Load R packages ##
#####################
required_libraries <- c("data.table", "tidyverse", "purrr", "optparse", "preprocessCore")
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
        c("-n", "--network_dir"),
        type = "character",
        default = NULL,
        help = "Path to the directory containing network RData files.",
        metavar = "character"
    ),
    optparse::make_option(
        c("-s", "--sample_file"),
        type = "character",
        default = NULL,
        help = "Path to the sample file containing cancer type information.",
        metavar = "character"
    ),
    optparse::make_option(
        c("-o", "--output_dir"),
        type = "character",
        default = NULL,
        help = "Path to the output directory for saving normalized network RData files.",
        metavar = "character"
    )
)

# Parse command line arguments
opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

# Assign parsed arguments to variables
NETWORK_DIR <- opt$network_dir
OUTPUT_DIR <- opt$output_dir
SAMPLE_FILE <- opt$sample_file

# Source functions
source("workflow/bin/utils_pcp_fn.R")

cat("Starting network quantile normalization...\n")
cat(sprintf("Network directory: %s\n", NETWORK_DIR))
cat(sprintf("Sample file: %s\n", SAMPLE_FILE))
cat(sprintf("Output directory: %s\n", OUTPUT_DIR))


####################
## Load data      ##
####################

# Ensure output directory exists
if (!dir.exists(OUTPUT_DIR)) {
    cat("Creating output directory:", OUTPUT_DIR, "\n")
    dir.create(OUTPUT_DIR, recursive = TRUE)
}


# Load sample file
if (!file.exists(SAMPLE_FILE)) {
    stop(sprintf("Cannot find sample file: %s", SAMPLE_FILE))
}
samples <- fread(SAMPLE_FILE)

# Get unique cancer types
cancers <- unique(samples$cancer)
cancers <- tolower(gsub("TCGA-", "", cancers))
cat(sprintf("Found %d cancer types: %s\n", length(cancers), paste(head(cancers), collapse = ", ")))

####################
## Process networks ##
####################

cat("Processing networks by cancer type...\n")
for (i in 1:length(cancers)) {
    cancer <- cancers[i]
    cat(sprintf("\nProcessing cancer type %d/%d: %s\n", i, length(cancers), cancer))
    
    # Quantile normalize networks
    tryCatch({
        net_norm <- quantile_normalize_net(cancer, NETWORK_DIR)
        output_file <- file.path(OUTPUT_DIR, paste0("net_norm_TCGA-", toupper(cancer), ".RData"))
        cat(sprintf("  Saving normalized network: %s\n", output_file))
        save(net_norm, file = output_file)
        cat(sprintf("  Successfully processed %s\n", cancer))
        rm(net_norm)
        gc()
    }, error = function(e) {
        cat(sprintf("  Error processing %s: %s\n", cancer, e$message))
    })
}
