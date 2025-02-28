required_libraries <- c("data.table", "dplyr", "doParallel",
            "optparse", "survminer", "survival", 
            "magrittr", "gridExtra", "purrr")

for (lib in required_libraries) {
  suppressPackageStartupMessages(library(lib, character.only = TRUE,
                                quietly = TRUE))
}

### Options
options(stringsAsFactors = FALSE)

### Command line options
option_list <- list(
    optparse::make_option(
        c("-t", "--tumor"),
        type = "character",
        default = NULL,
        help = "cancer type.",
        metavar = "character"),
    optparse::make_option(
        c("-c", "--clinical"),
        type = "character",
        default = NULL,
        help = "Path to the clinical file.",
        metavar = "character"),
    optparse::make_option(
        c("-i", "--input_dir"),
        type = "character",
        default = NULL,
        help = "Path to the pd1_dir directory.",
        metavar = "character"),
    optparse::make_option(
        c("-o", "--output_dir"),
        type = "character",
        default = NULL,
        help = "Path to the output_dir directory.",
        metavar = "character")
)

opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

## Initialize variable
INPUT_TUMOR <- opt$tumor
CLINICAL_FILE <- opt$clinical
INPUT_DIR <- opt$input_dir
OUTPUT_DIR <- opt$output_dir

## Debug
# CANCER <- c("acc")
# CLINICAL_FILE  <- "data/clinical_curated/NIHMS978596-supplement-1.xlsx"
# INPUT_DIR <- "data/pd1_data/"
# OUTPUT_DIR <- "intermediate_results/cox_results/"

# Debugging: Print parsed arguments
message("Parsed arguments:")
message("Tumor: ", INPUT_TUMOR)
message("Clinical File: ", CLINICAL_FILE)
message("Input Directory: ", INPUT_DIR)
message("Output Directory: ", OUTPUT_DIR)


## Functions
source("bin/cox_regression_fn.R")
## Create output directory
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
clin <- load_clin_curated(CLINICAL_FILE)
res <- run_glmnet_ntimes_pd1(tumor = INPUT_TUMOR, 
                            clin = clin,
                            pd1_dir = INPUT_DIR,
                            ind_scores_dir = NULL, 
                            pathways = NULL,
                            number_folds = 5,
                            ntimes = 1,
                            ncores = 1,
                            alpha = 0.3)

write.table(res,
                file.path(OUTPUT_DIR, paste0(INPUT_TUMOR, "_cox_results.txt")),
                col.names = T, row.names =  F, sep = "\t", quote = F)

