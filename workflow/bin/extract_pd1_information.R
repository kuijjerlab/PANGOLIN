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
        help = "Cancer type to analyze (e.g., BRCA, LUAD, UVM).",
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
        help = "Path to pathway gene set file (GMT format).",
        metavar = "character"
    ),
    optparse::make_option(
        c("-z", "--pd1_edges_file"),
        type = "character",
        default = NULL,
        help = "Path to PD-1 edges file.",
        metavar = "character"
    ),
    optparse::make_option(
        c("-r", "--pd1_net_file"),
        type = "character",
        default = NULL,
        help = "Path to PD-1 network file.",
        metavar = "character"
    ))


# Parse command line arguments
opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)
TUMOR_TYPE <- opt$tumor_type
NETWORK_FILE <- opt$network_file
EDGE_FILE <- opt$edge_file
PATHWAY_FILE <- opt$pathway_file

PD1_EDGES_FILE <- opt$pd1_edges_file
PD1_NET_FILE <- opt$pd1_net_file
cat("Starting PD-1 edge and network extraction...\n")
cat(sprintf("Tumor type: %s\n", TUMOR_TYPE))
cat(sprintf("Network file: %s\n", NETWORK_FILE))
cat(sprintf("Edge file: %s\n", EDGE_FILE))
cat(sprintf("Pathway file: %s\n", PATHWAY_FILE))
cat(sprintf("PD-1 edges output file: %s\n", PD1_EDGES_FILE))
cat(sprintf("PD-1 network output file: %s\n", PD1_NET_FILE))
####################
## Load data      ##
####################
source("workflow/bin/utils_pcp_fn.R")
cat("Extracting PD-1 related edges and network...\n")

extract_pd1_edges_norm(network_files = NETWORK_FILE,
                        cancer = TUMOR_TYPE,
                        edge_file = EDGE_FILE,
                        pathway_file = PATHWAY_FILE,
                        pd1_edges_file = PD1_EDGES_FILE,
                        pd1_net_file = PD1_NET_FILE)