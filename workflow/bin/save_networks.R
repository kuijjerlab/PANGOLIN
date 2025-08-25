#####################
## Load R packages ##
#####################
required_libraries <- c("gtools", "purrr", "tidyverse", "data.table", "optparse", "plyr")
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
        help = "Path to the directory containing LIONESS network txt files.",
        metavar = "character"
    ),
    optparse::make_option(
        c("-i", "--lioness_sample_mapping"),
        type = "character",
        default = NULL,
        help = "Path to the LIONESS sample mapping file.",
        metavar = "character"
    ),
    optparse::make_option(
        c("-o", "--output_dir"),
        type = "character",
        default = NULL,
        help = "Path to the output directory for saving combined network RData files.",
        metavar = "character"
    )
)

# Parse command line arguments
opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

# Assign parsed arguments to variables
NETWORK_DIR <- opt$network_dir
LIONESS_SAMPLE_MAPPING <- opt$lioness_sample_mapping
OUTPUT_DIR <- opt$output_dir

cat("Starting LIONESS network combination and saving...\n")
cat(sprintf("Network directory: %s\n", NETWORK_DIR))
cat(sprintf("Sample mapping file: %s\n", LIONESS_SAMPLE_MAPPING))
cat(sprintf("Output directory: %s\n", OUTPUT_DIR))

####################
## Load data      ##
####################

# Source helper functions
source("workflow/bin/split_save_networks_fn.R")

# Load LIONESS network manifest
cat("Loading LIONESS sample mapping...\n")
info_net <- fread(LIONESS_MAPPING_FILE)

# Get unique cancer types
cancers <- unique(info_net$cancer)

####################
## Process networks ##
####################

# Change to network directory for relative path processing
setwd(NETWORK_DIR)

# Process each cancer type
cat("Processing networks by cancer type...\n")
for (i in 1:length(cancers)) {
    cancer_type <- cancers[i]
    cat(sprintf("\nProcessing cancer type %d/%d: %s\n", i, length(cancers), cancer_type))
    # Get networks for this cancer type
    cancer_networks <- info_net[info_net$cancer == cancer_type, ]
    cat(sprintf("  Found %d networks for %s\n", nrow(cancer_networks), cancer_type))
    # Combine and save networks
    tryCatch({
        save_combined_networks(cancer_type, info_net, OUTPUT_DIR)
        cat(sprintf("  Successfully processed %s\n", cancer_type))
    }, error = function(e) {
        cat(sprintf("  Error processing %s: %s\n", cancer_type, e$message))
    })
}


