#####################
## Load R packages ##
#####################
# List of packages to be loaded
required_libraries <- c("data.table", "optparse", "tidyr", "dplyr",
                "cola", "stringr", "tidyverse", "viridis", "ggplot2")

for (lib in required_libraries) {
  suppressPackageStartupMessages(library(lib, character.only = TRUE,
                                quietly = TRUE))
}
set.seed(1234)
####################
## Read arguments ##
####################
option_list = list(
    make_option(c("-t", "--tumor"), type = "character", default = NULL,
            help = "Tumor type to filter (e.g., ACC, BRCA, LUNG)",
            metavar = "character"),
    make_option(c("-e", "--exp_file"), type = "character", default = NULL,
            help = "Path to the expression file",
            metavar = "character"),
    make_option(c("-s", "--samples_file"), type = "character", default = NULL,
            help = "Path to the samples file",
            metavar = "character"),
    make_option(c("-i", "--indegree_dir"), type = "character", default = NULL,
            help = "Path to the indegree directory for a specific cancer",
            metavar = "character"),
    make_option(c("-d", "--datatype"), type = "character", default = NULL,
            help = "Either indegree or expression",
            metavar = "character"),
    make_option(c("-n", "--number_cores"), type = "integer", default = NULL,
            help = "Number of cores to use. Default is 1",
            metavar = "integer"),
    make_option(c("-v", "--top_value_method"), type = "character", 
            default = NULL,
            help = "Top value method in cola clustering.
            Default is ATC. Other options check in cola package",
            metavar = "character"),
    make_option(c("-p", "--partition_method"), type = "character",
            default = NULL,
            help = "Partition method in cola clustering,
            Default is kmeans. Other options check in cola package",
            metavar = "character"),
    make_option(c("-k", "--max_k"), type = "integer", default = NULL,
            help = "Maximum number of clusters to test. Default is 6",
            metavar = "integer"),
    make_option(c("-o", "--output"), type = "character", default = NULL,
            help = "Output directory",
            metavar = "character"))


opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)



## Initialize variable
TUMOR <- opt$tumor
EXPRESSION_FILE <- opt$exp_file
SAMPLES_FILE <- opt$samples_file
IND_DIR <- opt$indegree_dir
DATATYPE <- opt$datatype
N_CORES <- opt$number_cores
TOP_VALUE_METHOD <- opt$top_value_method
PARTITION_METHOD <- opt$partition_method
MAX_K <- opt$max_k
OUTPUT_DIR <- opt$output

source("bin/cola_clustering_fn.R")
# perform conssensus clustering on indegree and save the results
perform_cola_clustering(cancer = TUMOR,
                                samples_file = SAMPLES_FILE, 
                                exp_file = EXPRESSION_FILE,
                                indegree_dir = IND_DIR,
                                datatype = DATATYPE,
                                n_cores = N_CORES,
                                top_value_method = TOP_VALUE_METHOD,
                                partition_method = PARTITION_METHOD,
                                max_k = MAX_K,
                                p_sampling = 0.8,
                                partition_repeat = 1000,
                                output_dir = OUTPUT_DIR
                            )
