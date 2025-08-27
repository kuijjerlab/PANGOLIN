#####################
## Load R packages ##
#####################
required_libraries <- c("data.table", "optparse", "tidyverse", "gtools", "purrr")
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
        c("-t", "--tumor_type"),
        type = "character",
        default = NULL,
        help = "Cancer type to analyze (e.g., BRCA, LUAD, UVM).",
        metavar = "character"
    ),
    optparse::make_option(
        c("-n", "--network_dir"),
        type = "character",
        default = NULL,
        help = "Path to directory containing quantile-normalized network RData files (e.g., net_norm_TCGA-BRCA.RData).",
        metavar = "character"
    ),
    optparse::make_option(
        c("-e", "--edge_file"),
        type = "character",
        default = NULL,
        help = "Path to edge list file defining network structure (TF-gene regulatory relationships).",
        metavar = "character"
    ),
    optparse::make_option(
        c("-o", "--output_file"),
        type = "character",
        default = NULL,
        help = "Path for output RData file containing calculated indegree matrix for the specified cancer type.",
        metavar = "character"
    )
)

# Parse command line arguments
opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

# Assign parsed arguments to variables
TUMOR_TYPE <- opt$tumor_type
NETWORK_DIR <- opt$network_dir
EDGE_FILE <- opt$edge_file
OUTPUT_FILE <- opt$output_file

cat("Starting indegree calculation...\n")
cat(sprintf("Cancer type: %s\n", TUMOR_TYPE))
cat(sprintf("Network directory: %s\n", NETWORK_DIR))
cat(sprintf("Edge file: %s\n", EDGE_FILE))
cat(sprintf("Output file: %s\n", OUTPUT_FILE))

####################
## Load functions ##
####################
# Source utility functions for indegree calculation
source("workflow/bin/utils_pcp_fn.R")

####################
## Calculate indegree ##
####################
cat("Calculating normalized indegree for cancer type:", TUMOR_TYPE, "\n")

# Calculate indegree using normalized networks
# This function loads the normalized network for the specified cancer type
# and calculates the indegree (number of incoming edges) for each gene
ind <- calculate_indegree_norm(TUMOR_TYPE, NETWORK_DIR, EDGE_FILE)

cat("Saving indegree results to:", OUTPUT_FILE, "\n")
save(ind, file = OUTPUT_FILE)
cat("Indegree calculation completed successfully.\n")
