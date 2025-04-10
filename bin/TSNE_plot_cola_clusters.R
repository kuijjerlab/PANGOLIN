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
        metavar = "character"),
    make_option(
        c("-o", "--output_dir"),
        type = "character",
        default = NULL,
        help = "Path to the output directory.",
        metavar = "character"),
    make_option(
        c("-f", "--output_figure_dir"),
        type = "character",
        default = NULL,
        help = "Path to the output figure directory.",
        metavar = "character"))


set.seed(1234)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

## Initialize variable
TUMOR_DIR_MAIN <- opt$tumor_dir
OUTPUT_DIR <- opt$output_dir
BEST_K_IND <- opt$best_cola_k_indegree
BEST_K_EXP <- opt$best_cola_k_expression
OUTPUT_DIR <- opt$output_dir
OUTPUT_FIG_DIR <- opt$output_figure_dir

source("bin/cola_clustering_fn.R")

res_files_indegree <- 
            list.files(TUMOR_DIR_MAIN, 
            full = TRUE, 
            recursive = TRUE, 
            pattern = "results_indegree")

# extract the clucter information for the best number of clusters 
# plot TSNE plots for all cancers
res_ind <- extract_relevant_classes(res_files_indgree, BEST_K_IND)
tumors <- unique(res_ind$cancer)
tsne_res_all_ind <- perform_tsne_all_cancers(tumors, res_files_indegree)
write.table(tsne_res_all_ind, 
            file.path(OUTPUT_DIR, "tsne_ind_cancers.txt"),
            col.names = T, row.names = F, sep = "\t", quote = F)

res_ind$Dim1 <- tsne_res_all_ind$Dim1[match(res_ind$ID, tsne_res_all_ind$id)]
res_ind$Dim2 <- tsne_res_all_ind$Dim2[match(res_ind$ID, tsne_res_all_ind$id)]
res_ind$class <- paste("cl", res_ind$class, sep = "_")

pdf(file.path(OUTPUT_FIG_DIR, "cola_clusters_indegree.pdf"), 
            width = 10, height = 14)
p1 <- plot_tsne_clusters_cola(res_ind)
print(p1)
dev.off()

# EXPRESSION

res_files_expression <- 
            list.files(TUMOR_DIR_MAIN, 
            full = TRUE, 
            recursive = TRUE, 
            pattern = "results_expression")


res_exp <- extract_relevant_classes(res_files_expression, BEST_K_EXP)
tumors <- unique(res_exp$cancer)
tsne_res_all_exp <- perform_tsne_all_cancers(tumors, res_files_expression)
write.table(tsne_res_all_exp, 
            file.path(OUTPUT_DIR, "tsne_exp_cancers.txt"),
            col.names = T, row.names = F, sep = "\t", quote = F)

res_exp$Dim1 <- tsne_res_all_exp$Dim1[match(res_exp$ID, tsne_res_all_exp$id)]
res_exp$Dim2 <- tsne_res_all_exp$Dim2[match(res_exp$ID, tsne_res_all_exp$id)]
res_exp$class <- paste("cl", res_exp$class, sep = "_")

pdf(file.path(OUTPUT_FIG_DIR, "cola_clusters_expression.pdf"),
            width = 10, height = 14)
p1 <- plot_tsne_clusters_cola(res_exp)
print(p1)
dev.off()
