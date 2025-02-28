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
        c("-tf", "--tumor_file"),
        type = "character",
        default = NULL,
        help = "Path to the file containing all tumors",
        metavar = "character"),
    optparse::make_option(
        c("-c", "--clinical"),
        type = "character",
        default = NULL,
        help = "Path to the clinical file.",
        metavar = "character"),
    optparse::make_option(
        c("-o", "--output_dir"),
        type = "character",
        default = NULL,
        help = "Path to the output_dir directory.",
        metavar = "character"))

opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

## Initialize variable
TUMOR_FILE <- pt$tumor_file
CLINICAL_FILE <- opt$clinical
OUTPUT_DIR <- opt$output_dir
