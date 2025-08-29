
#####################
## Load R packages ##
#####################
# List of packages to be loaded
required_libraries <- c("TCGAbiolinks", "SummarizedExperiment", 
            "dplyr", "DT", "optparse")

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
    make_option(c("-t", "--tumor"), type = "character", default = NULL,
            help = "Tumor type to filter (e.g., ACC, BRCA, LUNG)",
            metavar = "character"),
    make_option(c("-o", "--output_file"), type = "character", default = NULL,
            help = "Path to the output file",
            metavar = "character")
)

opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

## Initialize variable ##

TUMOR_TYPE <- opt$tumor
OUTPUT_FILE <- opt$output_file

source("workflow/bin/download_gdc_fn.R")

# Download expression data for the specified tumor type #
project <- paste0("TCGA-", TUMOR_TYPE)
download_gdc_expression(project, OUTPUT_FILE)

