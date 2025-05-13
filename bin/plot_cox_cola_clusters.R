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
    make_option(c("-r", "--cox_results_cluster_file"), type = "character", 
            default = NULL,
            help = "Path to the file with the Cox results for cola clusters",
            metavar = "character"),
    make_option(c("-c", "--cancer_color_file"),
            type = "character", default = NULL,
            help = "Path to the cancer color file",
            metavar = "character"),
    make_option(c("-o", "--output_file"), type = "character", 
            default = NULL,
            help = "Path to the output figure file",
            metavar = "character"))



# Parse the command-line arguments
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Assign parsed arguments to variables
COX_RESULTS_COLA_FILE <- opt$cox_results_cluster_file
OUTPUT_FIGURE_FILE <- opt$output_file
CANCER_COLOR_FILE <- opt$cancer_color_file

########################
## Load Helper Scripts ##
########################
# source required functions
source("bin/cox_regression_tumor_fn.R")

######################
## Load Data #########
######################
df <- fread(COX_RESULTS_COLA_FILE)
pfi_cancer <-  c("brca", "lgg", "prad", "read", "tgct", "thca", "thym")
pfi_cancer <- toupper(pfi_cancer)
# Filter the results based on the type and cancer
df <- df %>%
            filter(type == "PFI" & cancer %in% pfi_cancer |
                    type == "OS" & !cancer %in% pfi_cancer)

df$log10pval <- -log10(df$pval)
df$log10pval <- -log10(df$pval)
cancer_colors <- fread(CANCER_COLOR_FILE)
cancer_colors$cancer_type <- toupper(cancer_colors$cancer_type)
df$cancer <- toupper(df$cancer)
df$colour <- cancer_colors$colour[match(df$cancer, cancer_colors$cancer_type)]
df$shape <- cancer_colors$shape[match(df$cancer, cancer_colors$cancer_type)]
df$cluster <- df$k
plot_cox_results_cola_clusters(df, cancer_colors)

########################################
## Plot Cox Results for Cola Clusters ##
########################################

pdf(OUTPUT_FIGURE_FILE, width = 8, height = 8)
plot_cox_results_cola_clusters(df, cancer_colors)
dev.off()
