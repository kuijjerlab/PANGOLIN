#####################
## Load R packages ##
#####################
required_libraries <- c("data.table", "optparse", "tidyr", "dplyr",
                        "cola", "stringr", "tidyverse", "viridis", "ggplot2")

for (lib in required_libraries) {
  suppressPackageStartupMessages(library(lib, 
    character.only = TRUE, quietly = TRUE))
}
set.seed(1234)

####################
## Read arguments ##
####################
option_list = list(
    make_option(c("--tumor"), type = "character", default = NULL,
                help = "Tumor type to filter (e.g., ACC, BRCA)",
                metavar = "character"),
    make_option(c("--exp_file"), type = "character", default = NULL,
                help = "Path to the expression file",
                metavar = "character"),
    make_option(c("--samples_file"), type = "character", default = NULL,
                help = "Path to the samples file",
                metavar = "character"),
    make_option(c("--indegree_file"), type = "character", default = NULL,
                help = "Path to the indegree file for a specific cancer",
                metavar = "character"),
    make_option(c("--datatype"), type = "character", default = NULL,
                help = "Either indegree or expression",
                metavar = "character"),
    make_option(c("--number_cores"), type = "integer", default = NULL,
                help = "Number of cores to use. Default is 1",
                metavar = "integer"),
    make_option(c("--top_value_method"), type = "character", 
                default = NULL,
                help = "Top value method in cola clustering. Default is ATC.",
                metavar = "character"),
    make_option(c("--partition_method"), type = "character",
                default = NULL,
                help = "Partition method in cola clustering.
                Default is kmeans.",
                metavar = "character"),
    make_option(c("--max_k"), type = "integer", default = NULL,
                help = "Maximum number of clusters to test. Default is 6",
                metavar = "integer"),
    make_option(c("--output_best_k"), type = "character", default = NULL,
                help = "Path to save the best k information",
                metavar = "character"),
    make_option(c("--output_results"), type = "character", default = NULL,
                help = "Path to save the main results object",
                metavar = "character"),
    make_option(c("--output_membership"), type = "character", default = NULL,
                help = "Path to save the membership information",
                metavar = "character"),
    make_option(c("--output_statistics"), type = "character", default = NULL,
                help = "Path to save the statistics",
                metavar = "character"),
    make_option(c("--output_classes"), type = "character", default = NULL,
                help = "Path to save the classes information",
                metavar = "character"),
    make_option(c("--output_collected_plots"), type = "character", 
                default = NULL,
                help = "Path to save the collected plots PDF",
                metavar = "character"),
    make_option(c("--output_tsne_pdf"), type = "character", default = NULL,
                help = "Path to save the t-SNE plot PDF",
                metavar = "character"),
    make_option(c("--output_partition_pdf"), type = "character", default = NULL,
                help = "Path to save the partition selection plot PDF",
                metavar = "character")
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

## Initialize variables
TUMOR <- opt$tumor
EXPRESSION_FILE <- opt$exp_file
SAMPLES_FILE <- opt$samples_file
INDEGREE_FILE <- opt$indegree_file
DATATYPE <- opt$datatype
N_CORES <- opt$number_cores
TOP_VALUE_METHOD <- opt$top_value_method
PARTITION_METHOD <- opt$partition_method
MAX_K <- opt$max_k

OUTPUT_BEST_K <- opt$output_best_k
OUTPUT_RESULTS <- opt$output_results
OUTPUT_MEMBERSHIP <- opt$output_membership
OUTPUT_STATISTICS <- opt$output_statistics
OUTPUT_CLASSES <- opt$output_classes
OUTPUT_COLLECTED_PLOTS <- opt$output_collected_plots
OUTPUT_TSNE_PDF <- opt$output_tsne_pdf
OUTPUT_PARTITION_PDF <- opt$output_partition_pdf

source("workflow/bin/cola_clustering_fn.R")

res <- perform_cola_clustering(cancer = TUMOR,
                                samples_file = SAMPLES_FILE,
                                exp_file = EXPRESSION_FILE,
                                indegree_file = INDEGREE_FILE,
                                datatype = DATATYPE,
                                n_cores = N_CORES,
                                top_value_method = TOP_VALUE_METHOD,
                                partition_method = PARTITION_METHOD,
                                max_k = MAX_K,
                                p_sampling = 0.8,
                                partition_repeat = 1000
                            )

save_results(res = res, 
            cancer = TUMOR,
            datatype = DATATYPE,
            output_best_k =OUTPUT_BEST_K,
            output_results = OUTPUT_RESULTS,
            output_membership = OUTPUT_MEMBERSHIP,
            output_statistics = OUTPUT_STATISTICS,
            output_classes = OUTPUT_CLASSES,
            output_collected_plots = OUTPUT_COLLECTED_PLOTS,
            output_tsne_pdf = OUTPUT_TSNE_PDF,
            output_partition_pdf = OUTPUT_PARTITION_PDF)