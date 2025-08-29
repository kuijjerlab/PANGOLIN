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
    make_option(c("-k", "--clinical_file_tumor"), type = "character", 
                default = NULL,
                help = "Path to the clinical file for the tumor",
                metavar = "character"),
    make_option(c("-t", "--tumor_pd1_directory"), type = "character", 
                default = NULL,
                help = "Path to the PD1 directory for the tumor",
                metavar = "character"),
    make_option(c("-c", "--cluster_file_expression"), type = "character", 
                default = NULL,
                help = "Path to the cluster file for expression data",
                metavar = "character"),
    make_option(c("-m", "--cluster_file_indegree"), type = "character",
                default = NULL,
                help = "Path to the cluster file for indegree data",
                metavar = "character"),
    make_option(c("-o", "--output_file"), type = "character", 
                default = NULL,
                help = "Path to the output file to save the results",
                metavar = "character")
)

# Parse the command-line arguments
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Assign parsed arguments to variables
CLUSTER_FILE_EXP <- opt$cluster_file_expression
CLUSTER_FILE_IND <- opt$cluster_file_indegree
TUMOR_CLIN_FILE <- opt$clinical_file_tumor
TUMOR_PD1_DIR <- opt$tumor_pd1_directory
OUTPUT_FILE <- opt$output_file

########################
## Load Helper Scripts ##
########################
# source required functions
source("workflow/bin/cox_regression_tumor_fn.R")
source("workflow/bin/extract_clinical_data_fn.R")

########################
## Define Covariates ##
########################
# Specify covariates to include in the Cox regression model
COVARIATES_TO_USE <- c("gender", "age_at_initial_pathologic_diagnosis")


##############################
## Run Cox Regression Models ##
##############################
# Run univariate Cox regression for expression clusters
results_exp <- run_univariate_coxph_model(
    tumor_clin_file_path = TUMOR_CLIN_FILE,
    tumor_pd1_dir = TUMOR_PD1_DIR,
    covariates = COVARIATES_TO_USE,
    cluster_file = CLUSTER_FILE_EXP,
    datatype = "clusters",
    type_stat = "overall_pval"
)

# Run univariate Cox regression for indegree clusters
results_ind <- run_univariate_coxph_model(
    tumor_clin_file_path = TUMOR_CLIN_FILE,
    tumor_pd1_dir = TUMOR_PD1_DIR,
    covariates = COVARIATES_TO_USE,
    cluster_file = CLUSTER_FILE_IND,
    datatype = "clusters",
    type_stat = "overall_pval"
)

######################################
## Extract and Combine Cox Results ##
######################################
# Extract results for expression clusters
cox_res_exp <- extract_coxph_results_clusters(results_exp)
cox_res_exp$datatype <- "expression"

# Extract results for indegree clusters
cox_res_ind <- extract_coxph_results_clusters(results_ind)
cox_res_ind$datatype <- "indegree"

# Combine results into a single data frame
cox_res <- rbind(cox_res_exp, cox_res_ind)

##########################
## Save Results to File ##
##########################
# Save the combined Cox regression results to the specified output file
fwrite(cox_res, 
  file = OUTPUT_FILE, 
  sep = "\t", 
  row.names = FALSE,
  quote = FALSE)

