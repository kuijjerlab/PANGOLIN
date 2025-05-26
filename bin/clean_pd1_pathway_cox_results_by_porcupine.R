#####################
## Load R packages ##
#####################
# List of packages to be loaded
required_libraries <- c("circlize", "data.table", 
            "dplyr", "plyr", "ComplexHeatmap", "viridis")

for (lib in required_libraries) {
  suppressPackageStartupMessages(library(lib, character.only = TRUE,
                                quietly = TRUE))
}

####################
## Read arguments ##
####################
### Options
options(stringsAsFactors = FALSE)
### Command line options
option_list <- list(
    optparse::make_option(
        c("-c", "--cox_summary_all_cancers"),
        type = "character",
        default = NULL,
        help = "Path to the cox summary file.",
        metavar = "character"),
    optparse::make_option(
        c("-p", "--porcupine_filtered_results"),
        type = "character",
        default = NULL,
        help = "Path to the the porcupine filtered cancer results.",
        metavar = "character"),
    optparse::make_option(
        c("-f", "--cox_summary_all_cancers_filtered"),
        type = "character",
        default = NULL,
        help = "Path to the cox summary file filtered by porcupine results.",
        metavar = "character")
        )
opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

## Initialize variable
COX_SUMMARY_ALL <- opt$cox_summary_all_cancers
PORCUPINE_RESULTS <- opt$porcupine_filtered_results
COX_SUMMARY_FILTERED <- opt$cox_summary_all_cancers_filtered

##########################
## Load functions ##
##########################
source("bin/merge_patient_data_fn.R")

res <- clean_cox_results(cox_results_file = COX_SUMMARY_ALL,
                        porcupine_results_file = PORCUPINE_RESULTS, 
                        pval_threshold = 0.05)
write.table(res, file = COX_SUMMARY_FILTERED, 
            sep = "\t", quote = FALSE, row.names = FALSE)
##########################
