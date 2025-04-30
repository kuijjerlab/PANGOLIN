#####################
## Load R packages ##
#####################
# List of packages to be loaded
required_libraries <- c("data.table", "optparse", "tidyr", "dplyr",
                "cola", "stringr", "tidyverse", "viridis", "ggplot2",
                "ggsankey", "survival", "survminer", "rstatix",
                "gridExtra", "grid", "knitr")

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
        c("-f", "--figure_sanky_plot"),
        type = "character",
        default = NULL,
        help = "Path to the output figure of the sanky plots.",
        metavar = "character"))


set.seed(1234)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

## Initialize variable
TUMOR_DIR_MAIN <- opt$tumor_dir
BEST_K_IND <- opt$best_cola_k_indegree
BEST_K_EXP <- opt$best_cola_k_expression
FUGURE_SANKY <- opt$figure_sanky_plot

source("bin/cola_clustering_fn.R")
source("bin/cox_regression_tumor_fn.R")
source("bin/merge_patient_data_fn.R")
source("bin/sanky_plots_fn.R")


res_files_indegree <- 
            list.files(TUMOR_DIR_MAIN,
            full = TRUE,
            recursive = TRUE, 
            pattern = "results_indegree")
res_files_expression <- 
            list.files(TUMOR_DIR_MAIN,
            full = TRUE,
            recursive = TRUE, 
            pattern = "results_expression")


# extract relavant classes
cl_ids_ind <- extract_relevant_classes(res_files_indegree, BEST_K_IND)
cl_ids_exp <- extract_relevant_classes(res_files_expression, BEST_K_EXP)

best_k_ind <- load_and_prepare_data(BEST_K_IND)
best_k_ind <- best_k_ind %>%
                    group_by(cancer) %>%
                    slice_max(possible_clusters, n = 1) %>%
                    mutate(k = paste("k", possible_clusters, sep = "_")) %>%
                    mutate(cancer_k = paste(cancer, k, sep = "_"))
     
best_k_exp <- load_and_prepare_data(BEST_K_EXP)
best_k_exp <- best_k_exp %>%
                    group_by(cancer) %>%
                    slice_max(possible_clusters, n = 1) %>%
                    mutate(k = paste("k", possible_clusters, sep = "_")) %>%
                    mutate(cancer_k = paste(cancer, k, sep = "_"))

cl_ids_ind <- cl_ids_ind %>%
                    mutate(cancer_k = paste(cancer, k, sep = "_"))

cl_ids_exp <- cl_ids_exp %>%
                    mutate(cancer_k = paste(cancer, k, sep = "_"))
cl_ids_ind <- cl_ids_ind[cl_ids_ind$cancer_k %in% best_k_ind$cancer_k,]
cl_ids_exp <- cl_ids_exp[cl_ids_exp$cancer_k %in% best_k_exp$cancer_k,]

cl_ids_ind$datatype <- "indegree"
cl_ids_exp$datatype <- "expression"

cl_ids_all <- rbind(cl_ids_ind, cl_ids_exp)
cl_ids_all$bcr_patient_barcode <- make_bcr_code(cl_ids_all$ID)
cl_ids_all$cluster_id <- cl_ids_all$class

datasets <- make_dataset_for_sanky(cl_ids_all)
r_index <- datasets$ri_index
do.call("rbind", r_index)
datasets_to_plot <- datasets$datasets
plots <- lapply(datasets_to_plot, function(x) sanky_plot_viridis(x, x$cancer))
plots_per_page <- 36

# Calculate the number of pages needed
num_pages <- ceiling(length(plots) / plots_per_page)

# Create a function to subset plots for each page
get_plots_for_page <- function(page_number) {
  start_index <- (page_number - 1) * plots_per_page + 1
  end_index <- min(length(plots), start_index + plots_per_page - 1)
  return(plots[start_index:end_index])
}

pdf(FUGURE_SANKY, width = 10, height =10)
for (page in 1:num_pages) {
  plots_for_page <- get_plots_for_page(page)
  grid.arrange(grobs = plots_for_page, ncol = 6 , nrow = 6)
}
dev.off()
