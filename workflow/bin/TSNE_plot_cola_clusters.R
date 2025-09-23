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
        c("-i", "--indegree_files"),
        type = "character",
        default = NULL,
        help = "Space-separated list of files 
          containing results of indegree cola clustering.",
        metavar = "character"),
    make_option(
        c("-e", "--expression_files"),
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
        metavar = "character"),
 make_option(
        c("-q", "--figure_TSNE_indegree"),
        type = "character",
        default = NULL,
        help = "Path to the output TSNE figure for indegree.",
        metavar = "character"),
    make_option(
        c("-n", "--figure_TSNE_expression"),
        type = "character",
        default = NULL,
        help = "Path to the output TSNE figure for expression.",
        metavar = "character")
)


opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

## Initialize variables
INDEGREE_FILES <- unlist(strsplit(opt$indegree_files, " "))
EXPRESSION_FILES <- unlist(strsplit(opt$expression_files, " "))
BEST_K_IND <- opt$best_cola_k_indegree
BEST_K_EXP <- opt$best_cola_k_expression
FIG_TSNE_INDEGREE <- opt$figure_TSNE_indegree
FIG_TSNE_EXPRESSION <- opt$figure_TSNE_expression

source("workflow/bin/cola_clustering_fn.R")

# ##############################
# ## Process INDEGREE Results ##
# ##############################

if (length(INDEGREE_FILES) > 0) {
  cat("Processing indegree results...\n")
  
  # Extract cluster information for the best number of clusters
  res_ind <- extract_relevant_classes(INDEGREE_FILES, BEST_K_IND)
  tumors <- unique(res_ind$cancer)
  
  # Perform TSNE for all cancers
  tsne_res_all_ind <- perform_tsne_all_cancers(tumors, INDEGREE_FILES)
  res_ind$Dim1 <- tsne_res_all_ind$Dim1[match(res_ind$ID, tsne_res_all_ind$id)]
  res_ind$Dim2 <- tsne_res_all_ind$Dim2[match(res_ind$ID, tsne_res_all_ind$id)]
  res_ind$class <- paste("cl", res_ind$class, sep = "_")
  
  # Generate TSNE plot
  pdf(FIG_TSNE_INDEGREE, width = 10, height = 14)
  p1 <- plot_tsne_clusters_cola(res_ind)
  print(p1)
  dev.off()
  cat("TSNE plot for indegree saved to:", FIG_TSNE_INDEGREE, "\n")
} else {
  cat("No indegree files provided.\n")
}

################################
## Process EXPRESSION Results ##
################################
if (length(EXPRESSION_FILES) > 0) {
  cat("Processing expression results...\n")
  
  # Extract cluster information for the best number of clusters
  res_exp <- extract_relevant_classes(EXPRESSION_FILES, BEST_K_EXP)
  tumors <- unique(res_exp$cancer)
  
  # Perform TSNE for all cancers
  tsne_res_all_exp <- perform_tsne_all_cancers(tumors, EXPRESSION_FILES)
  res_exp$Dim1 <- tsne_res_all_exp$Dim1[match(res_exp$ID, tsne_res_all_exp$id)]
  res_exp$Dim2 <- tsne_res_all_exp$Dim2[match(res_exp$ID, tsne_res_all_exp$id)]
  res_exp$class <- paste("cl", res_exp$class, sep = "_")
  
  # Generate TSNE plot
  pdf(FIG_TSNE_EXPRESSION, width = 10, height = 14)
  p1 <- plot_tsne_clusters_cola(res_exp)
  print(p1)
  dev.off()
  cat("TSNE plot for expression saved to:", FIG_TSNE_EXPRESSION, "\n")
} else {
  cat("No expression files provided.\n")
}