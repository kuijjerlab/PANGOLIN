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
source("workflow/bin/cox_regression_tumor_fn.R")
source("workflow/bin/PRAD_clusters_analysis_fn.R")
source("workflow/bin/cola_clustering_fn.R")

datatype <- "clusters"
covariates <- c("gender", "age_at_initial_pathologic_diagnosis")
type_outcome <- c("PFI")
cluster_id <- "k_4"

cluster_indegree <- fread(CLUSTER_INDEGREE)
cat("clusters avaiable through this COLA run:\n")
unique(cluster_indegree$k)
clusters <- unique(cluster_indegree$k)
### Make a survival plot for the clusters
cat("Making survival plot for clusters:\n",  clusters[1], "\n")
if ("k_4" %in% clusters) {
        CLUSTER_ID <- "k_4"
    } else {
        CLUSTER_ID <- clusters[1]
    }
# Create survival plot and save it
survival_plot <- plot_PRAD_cox_fit(tumor_clin_file_path = TUMOR_CLIN_FILE,
                            covariates = covariates,
                            cluster_file = CLUSTER_INDEGREE,
                            datatype = datatype,
                            type_outcome = type_outcome,
                            cluster_id = CLUSTER_ID)

# Save survival plot
ggsave(OUTPUT_SURV_PLOT, survival_plot$plot, width = 10, height = 8, device = "pdf")
cat("Saved survival plot to:", OUTPUT_SURV_PLOT, "\n")

### Make a plot with fgsea results comparing cl1 to the other clusters

res <- get_degs_drgs_PRAD(tumor_clin_file_path = TUMOR_CLIN_FILE,
                              covariates = covariates,
                              cluster_file= CLUSTER_INDEGREE,
                              datatype = datatype,
                              type_outcome = type_outcome,
                              indegree_file = INDEGREE_FILE,
                              exp_file = EXP_FILE ,
                              samples_file = SAMPLES_FILE,
                              cluster_id =  CLUSTER_ID)

indegree_res <- run_fgsea_PRAD(res$indegree, GMT_FILE)
plots <- plot_fgsea_results(indegree_res)

# Handle single plot or multiple pages - but always create the expected output file
if (is.list(plots)) {
    # Multiple pages - save as multi-page PDF to the expected filename
    pdf(OUTPUT_FGSEA_PLOT, width = 10, height = 12)
    for (plot in plots) {
        print(plot)
    }
    dev.off()
    cat("Saved", length(plots), "pages to:", OUTPUT_FGSEA_PLOT, "\n")
} else {
    # Single plot
    ggsave(OUTPUT_FGSEA_PLOT, plots, width = 10, height = 12, device = "pdf")
    cat("Saved single plot to:", OUTPUT_FGSEA_PLOT, "\n")
}

