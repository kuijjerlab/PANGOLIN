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
        c("-t", "--tumor_clin_file_path"),
        type = "character",
        default = NULL,
        help = "Path to the clinical file.",
        metavar = "character"),
    optparse::make_option(
        c("-p", "--tumor_pd1_dir"),
        type = "character",
        default = NULL,
        help = "Path to the pd1 tumor directory.",
        metavar = "character"),
    optparse::make_option(
        c("-c", "--covariates"),
        type = "character",
        default = NULL,
        help = "covariates to include in the model.",
        metavar = "character"),
    optparse::make_option(
        c("-s", "--type_stat"),
        type = "character",
        default = "full_stats",
        help = "Type of statistics to report (full_stats, overall_pval).",
        metavar = "character"),
    optparse::make_option(
        c("-m", "--cox_model_summary"),
        type = "character",
        default = NULL,
        help = "Path to the output file of cox model summary.",
        metavar = "character"),
    optparse::make_option(
        c("-k", "--cox_predicted_risk"),
        type = "character",
        default = NULL,
        help = "Path to the output file of cox predicted risk scores.",
        metavar = "character")
        )

opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

## Initialize variable
TUMOR_CLIN_FILE <- opt$tumor_clin_file_path
TUMOR_PD1_DIR <- opt$tumor_pd1_dir
COVARIATES_TO_USE <- opt$covariates
DATATYPE <- opt$datatype
TYPE_STAT <- opt$type_stat
OUTPUT_FILE_COX_SUMMARY <- opt$cox_model_summary
OUTPUT_FILE_COX_SCORES <- opt$cox_predicted_risk

source("bin/cox_regression_tumor_fn.R")
dir.create(dirname(OUTPUT_FILE_COX_SUMMARY), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(OUTPUT_FILE_COX_SCORES), recursive = TRUE, showWarnings = FALSE)

COVARIATES_TO_USE = NULL
DATATYPE = "pd1_scores"
TYPE_STAT = "full_stats"

pd1_res <- run_univariate_coxph_model(tumor_clin_file_path = TUMOR_CLIN_FILE,
                                tumor_pd1_dir = TUMOR_PD1_DIR,
                                covariates = COVARIATES_TO_USE,
                                datatype = DATATYPE,
                                type_stat = TYPE_STAT)
pc_names <- c("PC1", "PC2")
pd1_res <- process_pd1_univarite_cox_res(pd1_res, pc_names)
write.table(pd1_res$coxph_model_data, OUTPUT_FILE_COX_SUMMARY,
    col.names = TRUE, sep = "\t", row.names = FALSE, quote = FALSE)

write.table(pd1_res$predicted_risk_data, OUTPUT_FILE_COX_SCORES,
    col.names = TRUE, sep = "\t", row.names = FALSE, quote = FALSE)

