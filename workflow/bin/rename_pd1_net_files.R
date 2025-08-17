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
        c("-p", "--path_pd1_network_file"),
        type = "character",
        default = NULL,
        help = "Path to pd1 network file",
        metavar = "character"),
    optparse::make_option(
        c("-i", "--path_to_indegree_file"),
        type = "character",
        default = NULL,
        help = "Path to indegree file",
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
PATH_PD1_NETWORK_FILE <- opt$path_pd1_network_file
PATH_TO_INDEGREE_FILE <- opt$path_to_indegree_file
OUTPUT_FILE <- opt$output_file


source("bin/rename_pd1_net_files_fn.R")
dir.create(dirname(OUTPUT_FILE), recursive = TRUE, showWarnings = FALSE)
pd1_net <- rename_pd1_net_file(PATH_TO_INDEGREE_FILE, PATH_PD1_NETWORK_FILE)
save(pd1_net, file = OUTPUT_FILE)

