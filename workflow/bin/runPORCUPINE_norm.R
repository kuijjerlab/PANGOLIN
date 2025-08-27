#####################
## Load R packages ##
#####################
required_libraries <- c("data.table", "optparse", "PORCUPINE")
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
        help = "Cancer type to analyze (e.g., BRCA, LUAD, UVM). Must match TCGA cancer codes.",
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
        c("-p", "--pathway_file"),
        type = "character",
        default = NULL,
        help = "Path to pathway gene set file (GMT format) defining gene sets for PORCUPINE enrichment analysis.",
        metavar = "character"
    ),
    optparse::make_option(
        c("-c", "--ncores"),
        type = "numeric",
        default = 1,
        help = "Number of CPU cores to use for parallel processing. Default: 1",
        metavar = "numeric"
    ),
    optparse::make_option(
        c("-o", "--output_directory"),
        type = "character",
        default = NULL,
        help = "Path to output directory for PORCUPINE results.",
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
PATHWAY_FILE <- opt$pathway_file
NCORES <- opt$ncores
OUTPUT_DIRECTORY <- opt$output_directory

cat("Starting PORCUPINE pathway analysis...\n")
cat(sprintf("Cancer type: %s\n", TUMOR_TYPE))
cat(sprintf("Network directory: %s\n", NETWORK_DIR))
cat(sprintf("Edge file: %s\n", EDGE_FILE))
cat(sprintf("Pathway file: %s\n", PATHWAY_FILE))
cat(sprintf("Number of cores: %s\n", NCORES))
cat(sprintf("Output directory: %s\n", OUTPUT_DIRECTORY))

####################
## Load functions ##
####################
# Source utility functions for PORCUPINE analysis
source("workflow/bin/utils_pcp_fn.R")

####################
## Run PORCUPINE ##
####################
cat("Running PORCUPINE pathway analysis for cancer type:", TUMOR_TYPE, "\n")

# Run PORCUPINE analysis on normalized networks
# This function performs pathway enrichment analysis using quantile-normalized
# LIONESS networks to identify pathway-level regulatory changes
runPORCUPINE_norm(TUMOR_TYPE,
                  NETWORK_DIR, 
                  EDGE_FILE,
                  PATHWAY_FILE,
                  OUTPUT_DIRECTORY,
                  NCORES)

cat("PORCUPINE analysis completed successfully.\n")


