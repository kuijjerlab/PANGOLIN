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
        c("-n", "--network_file"),
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
        c("-o", "--porcupine_results"),
        type = "character",
        default = NULL,
        help = "Path to the PORCUPINE results.",
        metavar = "character"
    ),
    optparse::make_option(
        c("-r", "--pathways_results_random"),
        type = "character",
        default = NULL,
        help = "Path to the directory containing randomization results.",
        metavar = "character"
    ),
    
    optparse::make_option(
        c("-k", "--individual_scores"),
        type = "character",
        default = NULL,
        help = "Path to the individual scores file.",
        metavar = "character"
    ),
    
    optparse::make_option(
        c("-g", "--pathways_results"),
        type = "character",
        default = NULL,
        help = "Path to the directory containing pathways results.",
        metavar = "character"
    )
)



# Parse command line arguments
opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

# Assign parsed arguments to variables
TUMOR_TYPE <- opt$tumor_type
NETWORK_FILE <- opt$network_file
EDGE_FILE <- opt$edge_file
PATHWAY_FILE <- opt$pathway_file
NCORES <- opt$ncores
PATHWAYS_RESULTS_FILE <- opt$pathways_results
PATHWAYS_RESULTS_RANDOM_FILE <- opt$pathways_results_random
PORCUPINE_RESULTS_FILE <- opt$porcupine_results
INDIVIDUAL_SCORES_FILE <- opt$individual_scores

cat("Starting PORCUPINE pathway analysis...\n")
cat(sprintf("Cancer type: %s\n", TUMOR_TYPE))
cat(sprintf("Network file: %s\n", NETWORK_FILE))
cat(sprintf("Edge file: %s\n", EDGE_FILE))
cat(sprintf("Pathway file: %s\n", PATHWAY_FILE))
cat(sprintf("Number of cores: %s\n", NCORES))
cat(sprintf("Pathways results file: %s\n", PATHWAYS_RESULTS_FILE))
cat(sprintf("Pathways random results file: %s\n", PATHWAYS_RESULTS_RANDOM_FILE))
cat(sprintf("PORCUPINE results file: %s\n", PORCUPINE_RESULTS_FILE))
cat(sprintf("Individual scores file: %s\n", INDIVIDUAL_SCORES_FILE))

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

runPORCUPINE_norm(cancer = TUMOR_TYPE,
                  network_files = NETWORK_FILE, 
                  edge_file = EDGE_FILE,
                  pathway_file = PATHWAY_FILE,
                  ncores_to_use = NCORES,
                  pathways_results_file = PATHWAYS_RESULTS_FILE,
                  pathways_results_random_file = PATHWAYS_RESULTS_RANDOM_FILE,
                  porcupine_results_file = PORCUPINE_RESULTS_FILE,
                  individual_scores_file = INDIVIDUAL_SCORES_FILE,
                  minSize = 5,
                  maxSize = 7)
cat("PORCUPINE analysis completed successfully.\n")


