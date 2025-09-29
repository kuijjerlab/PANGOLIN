#####################
## Load R packages ##
#####################
required_libraries <- c("gtools", "purrr", "tidyverse", "data.table", 
                        "optparse", "plyr")
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
        c("-f", "--network_files"),
        type = "character",
        default = NULL,
        help = "Space separated list of LIONESS network TXT file paths.",
        metavar = "character"
    ),
    optparse::make_option(
        c("-c", "--cancer_type"),
        type = "character",
        default = NULL,
        help = "Cancer type to process (e.g., BRCA, LUAD, etc.).",
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
        help = "Path to output directory for combined network RData files.",
        metavar = "character"
    )
)

# Parse command line arguments
opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

# Assign parsed arguments to variables
NETWORK_FILES <- opt$network_files
CANCER_TYPE <- opt$cancer_type
LIONESS_SAMPLE_MAPPING <- opt$lioness_sample_mapping
OUTPUT_DIR <- opt$output_dir

NETWORK_FILES <- unlist(strsplit(NETWORK_FILES, " "))
cat("Starting LIONESS network combination and saving...\n")
cat(sprintf("Network files: %s\n", NETWORK_FILES))
cat(sprintf("Cancer type: %s\n", CANCER_TYPE))
cat(sprintf("Sample mapping file: %s\n", LIONESS_SAMPLE_MAPPING))
cat(sprintf("Output directory: %s\n", OUTPUT_DIR))

####################
## Load data      ##
####################

source("workflow/bin/split_save_networks_fn.R")

# Parse the comma-separated list of network files
cat("Parsing LIONESS network file list...\n")
filelist <- trimws(unlist(strsplit(NETWORK_FILES, ",")))

if (length(filelist) == 0) {
    stop("Error: No LIONESS network files provided")
}

# Check if files exist
missing_files <- filelist[!file.exists(filelist)]
if (length(missing_files) > 0) {
    stop(sprintf("Error: The following files do not exist:\n%s", 
                 paste(missing_files, collapse = "\n")))
}

cat(sprintf("Found %d LIONESS network files for cancer type: %s\n", 
            length(filelist), CANCER_TYPE))

# Load LIONESS network manifest
cat("Loading LIONESS sample mapping...\n")
if (!file.exists(LIONESS_SAMPLE_MAPPING)) {
    stop(sprintf("Cannot find LIONESS sample mapping file: %s", 
                 LIONESS_SAMPLE_MAPPING))
}
info_net <- fread(LIONESS_SAMPLE_MAPPING)

# Filter mapping for the specific cancer type and provided files
cat(sprintf("Filtering mapping for cancer type: %s\n", CANCER_TYPE))
cancer_info_net <- info_net[info_net$cancer == CANCER_TYPE, ]

if (nrow(cancer_info_net) == 0) {
    stop(sprintf("No samples found for cancer type: %s", CANCER_TYPE))
}

# Filter to only include the files we actually have
file_basenames <- basename(filelist)
cancer_info_net <- cancer_info_net[cancer_info_net$file %in% file_basenames, ]

cat(sprintf("Found %d matching samples in mapping file\n", 
            nrow(cancer_info_net)))

####################
## Process networks ##
####################

# Process the specific cancer type
cat(sprintf("Processing cancer type: %s\n", CANCER_TYPE))
tryCatch({
    save_combined_networks(tumor = CANCER_TYPE, info_net = info_net, 
                          output_dir = OUTPUT_DIR)
    cat(sprintf("Successfully processed %s\n", CANCER_TYPE))
}, error = function(e) {
    stop(sprintf("Error processing %s: %s", CANCER_TYPE, e$message))
})



