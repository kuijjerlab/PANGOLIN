#####################
## Load R Packages ##
#####################
# List of required libraries
required_libraries <- c("data.table", "optparse", "tidyr", "dplyr",
                        "stringr", "tidyverse", "viridis", "ggplot2", "ggrepel")

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
                help = "Path to the file with individual-cluster information
                indegree (for PRAD)",
                metavar = "character"),
    optparse::make_option(
                c("-o", "--output_figure_file"), 
                type = "character",
                default = NULL,
                help = "Path to the output figure file (for PRAD)",
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
OUTPUT_FILE <- opt$output_figure_file

########################
## Load Helper Scripts ##
########################
# source required functions
source("bin/cox_regression_tumor_fn.R")


# CLUSTER_INDEGREE <- "/storage/kuijjerarea/tatiana/PANGOLIN/data_individual_cancers/PRAD/final_clusters/final_clusters_indegree_PRAD.txt"
# TUMOR_CLIN_FILE <-"/storage/kuijjerarea/tatiana/PANGOLIN/data_individual_cancers/PRAD/clinical/curated_clinical_PRAD.txt"
# TUMOR_PD1_DIR <- "/storage/kuijjerarea/tatiana/PANGOLIN/data_individual_cancers/PRAD/pd1_data"

datatype <- "clusters"
covariates <- c("gender", "age_at_initial_pathologic_diagnosis")


pdf(OUTPUT_FILE, width = 8, height = 8)

plot_PRAD_cox_fit(TUMOR_CLIN_FILE,
                              TUMOR_PD1_DIR,
                              covariates,
                              CLUSTER_INDEGREE,
                              datatype = c("clusters"))

dev.off()
