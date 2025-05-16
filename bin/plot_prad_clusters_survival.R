#####################
## Load R Packages ##
#####################
# List of required libraries
required_libraries <- c("data.table", "optparse", "tidyr", "dplyr",
                        "stringr", "tidyverse", "viridis", "ggplot2", "ggrepel",
                        "tidytext", "limma")

# Load each library and suppress startup messages
for (lib in required_libraries) {
  suppressPackageStartupMessages(library(lib, character.only = TRUE,
                                quietly = TRUE))
}

# Set a random seed for reproducibility
set.seed(1234)

####################
## Parse Arguments ##
####################
# Define command-line options
option_list <- list(
    optparse::make_option(
        c("-t", "--prad_clin_file_path"),
        type = "character",
        default = NULL,
        help = "Path to the clinical file (PRAD).",
        metavar = "character"),
    optparse::make_option(
        c("-p", "--prad_pd1_dir"),
        type = "character",
        default = NULL,
        help = "Path to the pd1 tumor directory (PRAD).",
        metavar = "character"),
    optparse::make_option(
        c("-m", "--prad_cluster_file_indegree"), 
        type = "character",
        default = NULL,
        help = "Path to the file with 
        individual-cluster information indegree (for PRAD)",
        metavar = "character"),
    optparse::make_option(
        c("-i", "--indegree_file"), 
        type = "character",
        default = NULL,
        help = "Path to the indegree file (for PRAD)",
        metavar = "character"),
    optparse::make_option(
        c("-e", "--exp_file"), 
        type = "character", 
        default = NULL,
        help = "Path to the expression file",
        metavar = "character"),
    optparse::make_option(
        c("-s", "--samples_file"), 
        type = "character",
        default = NULL,
        help = "Path to the samples file",
        metavar = "character"),
    optparse::make_option(
        c("-g", "--gmt_file"),
        type = "character",
        default = NULL,
        help = "Path to the GMT file for fgsea",
        metavar = "character"),
    optparse::make_option(
        c("-o", "--output_survival_plot"), 
        type = "character",
        default = NULL,
        help = "Path to the output survival plot (for PRAD)",
        metavar = "character"),
    optparse::make_option(
        c("-l", "--output_fgsea_plot"), 
        type = "character",
        default = NULL,
        help = "Path to the output fgsea figure (for PRAD)",
        metavar = "character")
        )
# Parse the command-line arguments
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Assign parsed arguments to variables
## Initialize variable
TUMOR_CLIN_FILE <- opt$prad_clin_file_path
TUMOR_PD1_DIR <- opt$prad_pd1_dir
CLUSTER_INDEGREE <- opt$prad_cluster_file_indegree
INDEGREE_FILE <- opt$indegree_file
EXP_FILE <- opt$exp_file
SAMPLES_FILE <- opt$samples_file
GMT_FILE <- opt$gmt_file
OUTPUT_SURV_PLOT <- opt$output_survival_plot
OUTPUT_FGSEA_PLOT <- opt$output_fgsea_plot


########################
## Load Helper Scripts ##
########################
# source required functions
source("bin/cox_regression_tumor_fn.R")
source("bin/PRAD_clusters_analysis_fn.R")
source("bin/cola_clustering_fn.R")

datatype <- "clusters"
covariates <- c("gender", "age_at_initial_pathologic_diagnosis")
type_outcome <- c("PFI")
cluster_id <- "k_4"

pdf(OUTPUT_SURV_PLOT, width = 8, height = 8)
plot_PRAD_cox_fit(TUMOR_CLIN_FILE,
                              TUMOR_PD1_DIR,
                              covariates,
                              CLUSTER_INDEGREE,
                              datatype,
                              type_outcome,
                              cluster_id)
dev.off()

### Make a plot with fgsea results comparing cl1 to the other clusters

res <- get_degs_drgs_PRAD(TUMOR_CLIN_FILE,
                              TUMOR_PD1_DIR,
                              covariates,
                              CLUSTER_INDEGREE,
                              datatype,
                              type_outcome,
                              cluster_id = "k_4",
                              indegree_file = INDEGREE_FILE,
                              exp_file = EXP_FILE ,
                              samples_file = SAMPLES_FILE)

indegree_res <- run_fgsea_PRAD(res$indegree, GMT_FILE)
pdf(OUTPUT_FGSEA_PLOT, width = 8, height = 8)
plot_fgsea_results(indegree_res)
dev.off()
