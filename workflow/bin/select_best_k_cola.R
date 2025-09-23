#####################
## Load R packages ##
#####################
required_libraries <- c("data.table", "optparse", "tidyr", "dplyr",
                        "cola", "stringr", "tidyverse", "viridis", "ggplot2")
for (lib in required_libraries) {
    suppressPackageStartupMessages(
        library(lib, character.only = TRUE, quietly = TRUE)
    )
}
####################
## Read arguments ##
####################
option_list = list(
    make_option(
        c("-i", "--indegree_best_k_files"),
        type = "character",
        default = NULL,
        help = "Space-separated list of files 
          containing results of indegree cola clustering.",
        metavar = "character"),
    make_option(
        c("-e", "--expression_best_k_files"),
        type = "character",
        default = NULL,
        help = "Space-separated list of files 
          containing results of expression cola clustering.",
        metavar = "character"),
    make_option(
        c("--best_cola_k_indegree"),
        type = "character",
        default = NULL,
        help = "Path to a file to save the 
          combined best k information for indegree.",
        metavar = "character"),
    make_option(
        c("--best_cola_k_expression"),
        type = "character",
        default = NULL,
        help = "Path to a file to save the 
          combined best k information for expression.",
        metavar = "character")
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

## Initialize variables
INDEGREE_K_FILES <- unlist(strsplit(opt$indegree_best_k_files, " "))
EXPRESSION_K_FILES <- unlist(strsplit(opt$expression_best_k_files, " "))

BEST_K_IND <- opt$best_cola_k_indegree
BEST_K_EXP <- opt$best_cola_k_expression

source("workflow/bin/cola_clustering_fn.R")

print(INDEGREE_K_FILES)
print(EXPRESSION_K_FILES)
length(INDEGREE_K_FILES)
length(EXPRESSION_K_FILES)
##############################
## Process INDEGREE Results ##
##############################
if (length(INDEGREE_K_FILES) > 0) {
  cat("Processing indegree best k files...\n")
  res_indegree_list <- lapply(INDEGREE_K_FILES, function(file) {
    load_result(file, object_name = "res_k")
  })
  res_all_indegree <- do.call(rbind, res_indegree_list)
  combined_res_indegree <- combine_k_results(res_all_indegree)
  write.table(combined_res_indegree, 
              BEST_K_IND,
              col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
  cat("Indegree results saved to:", BEST_K_IND, "\n")
} else {
  cat("No indegree files provided.\n")
}

################################
## Process EXPRESSION Results ##
################################
if (length(EXPRESSION_K_FILES) > 0) {
  cat("Processing expression best k files...\n")
  res_expression_list <- lapply(EXPRESSION_K_FILES, function(file) {
    load_result(file, object_name = "res_k")
  })
  res_all_expression <- do.call(rbind, res_expression_list)
  combined_res_expression <- combine_k_results(res_all_expression)
  write.table(combined_res_expression,
              BEST_K_EXP,
              col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
  cat("Expression results saved to:", BEST_K_EXP, "\n")
} else {
  cat("No expression files provided.\n")
}


