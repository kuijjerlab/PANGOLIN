#####################
## Load R packages ##
#####################

required_libraries <- c("data.table", "optparse", "tidyr", "dplyr",
                "cola", "stringr", "tidyverse", "viridis", "ggplot2")
for (lib in required_libraries) {
  suppressPackageStartupMessages(library(lib, character.only = TRUE,
                                quietly = TRUE))
}
####################
## Read arguments ##
####################
option_list = list(
    make_option(
        c("-d", "--tumor_dir"),
        type = "character",
        default = NULL,
        help = "Path to the the main tumor directory.",
        metavar = "character"),
   make_option(
        c("-i", "--best_cola_k_indegree"),
        type = "character",
        default = NULL,
        help = "Path to a file containing the best information 
                about selected number of clusters for indegree.",
        metavar = "character"),
    make_option(
        c("-e", "--best_cola_k_expression"),
        type = "character",
        default = NULL,
        help = "Path to a file containing the best information 
                about selected number of clusters for expression.",
        metavar = "character"))


opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

## Initialize variable
TUMOR_DIR_MAIN <- opt$tumor_dir
BEST_K_IND <- opt$best_cola_k_indegree
BEST_K_EXP <- opt$best_cola_k_expression

source("bin/cola_clustering_fn.R")
##### Load in INDEGREE cluster results
result_files_indegree <-
        list.files(TUMOR_DIR_MAIN, 
        recursive = TRUE,
        full = TRUE, pattern = "best_k_indegree")
res_indegree_list <- lapply(result_files_indegree, function(file) {
    load_result(file, object_name = "res_k")
})
res_all_indegree <- do.call(rbind, res_indegree_list)
combined_res_indegree <- combine_k_results(res_all_indegree)
write.table(combined_res_indegree, 
            BEST_K_IND,
            col.names = T, row.names = F, sep = "\t", quote = F)



#### Load in EXPRESSION cluster results
result_files_expression <-
        list.files(TUMOR_DIR_MAIN, 
        recursive = TRUE,
        full = TRUE, pattern = "best_k_expression")
res_expression_list <- lapply(result_files_expression, function(file) {
    load_result(file, object_name = "res_k")
})
res_all_expression <- do.call(rbind, res_expression_list)
combined_res_expression <- combine_k_results(res_all_expression)
write.table(combined_res_expression, 
            BEST_K_EXP,
            col.names = T, row.names = F, sep = "\t", quote = F)
