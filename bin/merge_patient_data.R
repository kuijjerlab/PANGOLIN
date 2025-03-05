#####################
## Load R packages ##
#####################

required_libraries <- c("data.table", "optparse", "tidyr", "dplyr")
for (lib in required_libraries) {
  suppressPackageStartupMessages(library(lib, character.only = TRUE,
                                quietly = TRUE))
}
####################
## Read arguments ##
####################
option_list <- list(
    optparse::make_option(
        c("-p", "--tumor_pd1_dir"),
        type = "character",
        default = NULL,
        help = "Path to the pd1 tumor directory.",
        metavar = "character"),
    optparse::make_option(
        c("-r", "--risk_score_file"),
        type = "character",
        default = NULL,
        help = "Path to the risk score file.",
        metavar = "character"),
    optparse::make_option(
        c("-i", "--immune_file"),
        type = "character",
        default = "",
        help = "File path to the cibersort deconvoluted immune data.",
        metavar = "character"),
    optparse::make_option(
        c("-o", "--output_file"),
        type = "character",
        default = NULL,
        help = "Output file path.",
        metavar = "character")
        )

opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

## Initialize variable ##

TUMOR_PD1_DIR <- opt$tumor_pd1_dir
RISK_SCORE_FILE <- opt$risk_score_file
IMMUNE_FILE <- opt$immune_file
OUTPUT_FILE <- opt$output_file

## Load functions ##
source("bin/cox_regression_tumor_fn.R")
source("bin/merge_patient_data_fn.R")

dir.create(dirname(OUTPUT_FILE), recursive = TRUE, showWarnings = FALSE)
res <- merge_patient_data(tumor_pd1_dir = TUMOR_PD1_DIR,
                    risk_score_file = RISK_SCORE_FILE,
                    immune_file = IMMUNE_FILE)

write.table(res, file = OUTPUT_FILE, 
    sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)