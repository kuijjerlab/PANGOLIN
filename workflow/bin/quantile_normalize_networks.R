#####################
## Load R packages ##
#####################
required_libraries <- c("data.table", "tidyverse", "purrr", "optparse",
                        "preprocessCore")
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
        c("-n", "--network_files"),
        type = "character",
        default = NULL,
        help = "Space separated list of all network file paths.",
        metavar = "character"
    ),
    optparse::make_option(
        c("-m", "--mapping_file"),
        type = "character",
        default = NULL,
        help = "Path to TSV mapping file with columns: net_file, cancer.",
        metavar = "character"
    ),
    optparse::make_option(
        c("-c", "--cancer"),
        type = "character",
        default = NULL,
        help = "Cancer type (e.g., BRCA, LUAD, etc.).",
        metavar = "character"
    ),
    optparse::make_option(
        c("-o", "--output_file"),
        type = "character",
        default = NULL,
        help = "Path to output file.",
        metavar = "character"
    )
)

# Parse command line arguments
opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

# Assign parsed arguments to variables
NETWORK_FILES <- opt$network_files
MAPPING_FILE <- opt$mapping_file
TUMOR <- opt$cancer
OUTPUT_FILE <- opt$output_file

# Source functions
source("workflow/bin/utils_pcp_fn.R")
NETWORK_FILES <- unlist(strsplit(NETWORK_FILES, " "))
cat("Starting network quantile normalization...\n")
cat(sprintf("Network files: %s\n", NETWORK_FILES))
cat(sprintf("Mapping file: %s\n", MAPPING_FILE))
cat(sprintf("Cancer type: %s\n", TUMOR))
cat(sprintf("Output file: %s\n", OUTPUT_FILE))


####################
## Process data   ##
####################


# Load the mapping filecat("Loading network-cancer mapping file...\n")
cat("Loading network-cancer mapping file...\n")
if (!file.exists(MAPPING_FILE)) {
    stop(sprintf("Mapping file does not exist: %s", MAPPING_FILE))
}

mapping <- fread(MAPPING_FILE)

# Get files for this specific cancer type
TUMOR <- toupper(TUMOR)
cancer_files <- mapping[mapping$cancer == TUMOR, ]$net_file

if (length(cancer_files) == 0) {
    stop(sprintf("No network files found for cancer type: %s", TUMOR))
    
}

cat(sprintf("Mapping found %d files for cancer %s: %s\n",
            length(cancer_files), TUMOR,
            paste(head(cancer_files, 3), collapse = ", ")))

# Filter input files to only those needed for this cancer
cancer_file_paths <- character(0)
for (net_file in cancer_files) {
    # Find the full path from input files that matches this filename
    matching_files <- NETWORK_FILES[basename(NETWORK_FILES) == net_file]
    if (length(matching_files) > 0) {
        cancer_file_paths <- c(cancer_file_paths, matching_files[1])
    } else {
        cat(sprintf("Warning: File %s not found in input files\n",
                    net_file))
    }
}

if (length(cancer_file_paths) == 0) {
    stop(sprintf("No input files match the required files for cancer %s",
                 TUMOR))
}

# Check if filtered files exist
missing_files <- cancer_file_paths[!file.exists(cancer_file_paths)]
if (length(missing_files) > 0) {
    stop(sprintf("Error: The following files do not exist:\n%s",
                 paste(missing_files, collapse = "\n")))
}

cat(sprintf("Selected %d network files for cancer type %s\n",
            length(cancer_file_paths), TUMOR))

# Ensure output directory exists
output_dir <- dirname(OUTPUT_FILE)
if (!dir.exists(output_dir)) {
    cat("Creating output directory:", output_dir, "\n")
    dir.create(output_dir, recursive = TRUE)
}

####################
## Process networks ##
####################

cat(sprintf("Processing cancer type: %s\n", TUMOR))
# Quantile normalize networks
tryCatch({
        net_norm <- quantile_normalize_net(cancer_file_paths)
        cat(sprintf("  Saving normalized network: %s\n", OUTPUT_FILE))
        save(net_norm, file = OUTPUT_FILE)
        cat(sprintf("  Successfully processed %s\n", TUMOR))
        rm(net_norm)
        gc()
    }, error = function(e) {
        cat(sprintf("  Error processing %s: %s\n", TUMOR, e$message))
    })
