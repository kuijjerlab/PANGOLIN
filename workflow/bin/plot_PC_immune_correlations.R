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
        c("-d", "--tumor_dir"),
        type = "character",
        default = NULL,
        help = "Path to the the main tumor directory.",
        metavar = "character"),
    optparse::make_option(
        c("-o", "--output_file"),
        type = "character",
        default = NULL,
        help = "Path to the output png figure file.",
        metavar = "character")
        )
opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

## Initialize variable
COX_SUMMARY_ALL <- opt$cox_summary_all_cancers
TUMOR_DIR_MAIN <- opt$tumor_dir
OUTPUT_PNG_FILE <- opt$output_file

source("workflow/bin/merge_patient_data_fn.R")
source("workflow/bin/plotting_fn.R")

res <- generate_PC_immune_correlation_table(cox_results_file = COX_SUMMARY_ALL,
                                cancer_dir = TUMOR_DIR_MAIN, 
                                pval_threshold = 0.05)
png(OUTPUT_PNG_FILE, width = 2000, height = 2000, res = 200) 
plot_pc_immune_correlations(res)
dev.off()

