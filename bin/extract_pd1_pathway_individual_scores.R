#####################
## Load R packages ##
#####################
required_libraries <- c("data.table", "dplyr", "doParallel",
            "survminer", "survival", "magrittr", "gridExtra",
            "purrr", "optparse")
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
        c("-t", "--tumor_pathways_mapping_path"),
        type = "character",
        default = NULL,
        help = "Path to tumor_pathways_mapping_path",
        metavar = "character"),
    optparse::make_option(
        c("-o", "--output_file"),
        type = "character",
        default = NULL,
        help = "Path to the output_file",
        metavar = "character"))

opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

## Initialize variable
TUMOR_PATH_MAPPING_FILE <- opt$tumor_pathways_mapping_path
OUTPUT_FILE <- opt$output_file


source("bin/cox_regression_tumor_fn.R")
dir.create(dirname(OUTPUT_FILE), recursive = TRUE, showWarnings = FALSE)
# TUMOR_PATH_MAPPING_FILE <-"data_individual_cancers/ESCA/porcupine/individual_scores_ESCA.RData"
# OUTPUT_FILE <- "data_individual_cancers/ESCA/pd1_data/pd1_individual_scores_norm_ESCA.RData"
ind_scores <- extract_pd1_pathway_individual_scores(TUMOR_PATH_MAPPING_FILE)
save(ind_scores, file = OUTPUT_FILE)



