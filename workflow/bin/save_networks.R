
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
        c("-f", "--network_dir"),
        type = "character",
        default = NULL,
        help = "Path to the directory containing LIONESS network TXT files.",
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

NETWORK_DIR <- opt$network_dir
LIONESS_SAMPLE_MAPPING <- opt$lioness_sample_mapping
OUTPUT_DIR <- opt$output_dir

filelist <- list.files(NETWORK_DIR, pattern = "^lioness\\.[0-9]+\\.txt$", full.names = TRUE)
cat("Starting LIONESS network combination and saving...\n")
cat(sprintf("Sample mapping file: %s\n", LIONESS_SAMPLE_MAPPING))
cat(sprintf("Output directory: %s\n", OUTPUT_DIR))

####################
## Load data      ##
####################

source("workflow/bin/split_save_networks_fn.R")

# Parse the comma-separated list of network files
cat("Parsing LIONESS network file list...\n")

if (length(filelist) == 0) {
    stop("Error: No LIONESS network files provided")
}

# Check if files exist
missing_files <- filelist[!file.exists(filelist)]
if (length(missing_files) > 0) {
    stop(sprintf("Error: The following files do not exist:\n%s", 
                 paste(missing_files, collapse = "\n")))
}

cat(sprintf("Found %d LIONESS network files in directory\n", length(filelist)))

# Load LIONESS network manifest
cat("Loading LIONESS sample mapping...\n")
if (!file.exists(LIONESS_SAMPLE_MAPPING)) {
    stop(sprintf("Cannot find LIONESS sample mapping file: %s", 
                 LIONESS_SAMPLE_MAPPING))
}
info_net <- fread(LIONESS_SAMPLE_MAPPING)

# Get all unique cancer types
cancer_types <- unique(info_net$cancer)
cat(sprintf("Found %d unique cancer types: %s\n", length(cancer_types), paste(cancer_types, collapse=", ")))

file_basenames <- basename(filelist)

for (CANCER_TYPE in cancer_types) {
    cat(sprintf("\nProcessing cancer type: %s\n", CANCER_TYPE))
    cancer_info_net <- info_net[info_net$cancer == CANCER_TYPE, ]
    if (nrow(cancer_info_net) == 0) {
        cat(sprintf("No samples found for cancer type: %s\n", CANCER_TYPE))
        next
    }
    # Filter to only include the files we actually have
    cancer_info_net <- cancer_info_net[cancer_info_net$file %in% file_basenames, ]
    cat(sprintf("Found %d matching samples in mapping file for %s\n", nrow(cancer_info_net), CANCER_TYPE))
    tryCatch({
        save_combined_networks(tumor = CANCER_TYPE, info_net = info_net, 
                              output_dir = OUTPUT_DIR)
        cat(sprintf("Successfully processed %s\n", CANCER_TYPE))
    }, error = function(e) {
        cat(sprintf("Error processing %s: %s\n", CANCER_TYPE, e$message))
    })
}



