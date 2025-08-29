
#####################
## Load R Packages ##
#####################
# List of required libraries
required_libraries <- c("data.table", "dplyr", "optparse",
                        "ggplot2", "ggrepel", "cowplot", 
                        "rstatix", "smplot2", "ggpubr")
# Load each library and suppress startup messages
for (lib in required_libraries) {
  suppressPackageStartupMessages(library(lib, character.only = TRUE,
                                quietly = TRUE))
}


####################
## Parse Arguments ##
####################
# Define command-line options
option_list <- list(
        make_option(c("-c", "--cox_results_file"),
                type = "character", 
                default = NULL,
                help = "Path to the cox results for the PD1 pathway",
                metavar = "character"),
        make_option(c("-r", "--results_pd1_groups"),
                type = "character",
                default = NULL,
                help = "Path to the result file comparing clinical groups",
                metavar = "character"),
        make_option(c("-s", "--results_pd1_numeric"),
                type = "character",
                default = NULL,
                help = "Path to the result file with numerical features",
                metavar = "character"),
        make_option(c("-f", "--tumor_main_dir"),
                type = "character", 
                default = NULL,
                help = "Path to the main directory containing all tumor types",
                metavar = "character"), 
        make_option(c("-o", "--output_figure_file"),
                type = "character",
                default = NULL,
                help = "Path to the output figure file",
                metavar = "character")
)

# Parse the command-line arguments
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Assign parsed arguments to variables
COX_RESULTS_FILE <- opt$cox_results_file
RESULTS_CATEGORICAL_FILE <- opt$results_pd1_groups
RESULTS_NUMERIC_FILE <- opt$results_pd1_numeric
TUMOR_MAIN_DIR <- opt$tumor_main_dir
OUTPUT_FIGURE <- opt$output_figure_file

########################
## Load Helper Scripts ##
########################
# source required functions
source("workflow/bin/clinical_association_pd1_fn.R")
source("workflow/bin/cox_regression_tumor_fn.R")
source("workflow/bin/extract_clinical_data_fn.R")
source("workflow/bin/plot_clinical_associations_fn.R")
set.seed(1234)



# plot the differences in PC scores between the clinical groups 

res_cat <- process_categorical_results(
                        res_categorical_file = RESULTS_CATEGORICAL_FILE,
                        coxph_results_file = COX_RESULTS_FILE)
res_num <- process_numeric_results(
                        res_numeric_file = RESULTS_NUMERIC_FILE,
                        coxph_results_file = COX_RESULTS_FILE)

data_all <- combine_data(res_cat, res_num)

results_categorical_all <- fread(RESULTS_CATEGORICAL_FILE)

data_all <- data_all %>%
      filter(padjust <= 0.01) %>%
      filter(abs(effect_size_cor) >= 0.4)

res_pc1 <- data_all %>%
      filter(principal_component == "PC1") %>%
      mutate(feature_cancer = paste0(cancer, "_", clin_feature))


features_pc1 <- unique(res_pc1$feature_cancer)
res_pc2 <- data_all %>%
      filter(principal_component == "PC2") %>%
      mutate(feature_cancer = paste0(cancer, "_", clin_feature))


features_pc2 <- unique(res_pc2$feature_cancer)

plots_pc1 <- list()

dirs <- list.dirs(TUMOR_MAIN_DIR, recursive = TRUE, full.names = TRUE)
pd1_dirs <- dirs[grepl("pd1_data", basename(dirs))]
clin_dirs <- dirs[grepl("clinical", basename(dirs))]

for (i in 1:length(features_pc1)) {
      splits <- split_string_function(features_pc1[i])
      tumor <- toupper(splits$cancer)
      feature_to_plot <- splits$feature
      pd1_dir_tumor <- pd1_dirs[grep(tumor, pd1_dirs)]
      tumor_clin_dir <- clin_dirs[grep(tumor, clin_dirs)]
      tumor_clin_file <- 
            list.files(tumor_clin_dir, recursive = TRUE, full.names = TRUE)
      plot <- plot_clin_feature(tumor = tumor,
                              results_categorical = results_categorical_all,
                              feature_to_plot = feature_to_plot,
                              component = "PC1",
                              clin_cancer_file = tumor_clin_file,
                              pd1_dir = pd1_dir_tumor)
      plots_pc1[[i]] <- plot
}


plots_pc2 <- list()

for (i in 1:length(features_pc2)) {
      splits <- split_string_function(features_pc2[i])
      tumor <- toupper(splits$cancer)
      feature_to_plot <- splits$feature
      pd1_dir_tumor <- pd1_dirs[grep(tumor, pd1_dirs)]
      tumor_clin_dir <- clin_dirs[grep(tumor, clin_dirs)]
      tumor_clin_file <- 
            list.files(tumor_clin_dir, recursive = TRUE, full.names = TRUE)
      plot <- plot_clin_feature(tumor = tumor,
                              results_categorical = results_categorical_all,
                              feature_to_plot = feature_to_plot,
                              component = "PC2",
                              clin_cancer_file = tumor_clin_file,
                              pd1_dir = pd1_dir_tumor)
      plots_pc2[[i]] <- plot
}

plots_all <- c(plots_pc1, plots_pc2)

# Determine the number of plots per page and calculate number of pages needed
plot_list <- plots_all
plots_per_page <- 9 
num_pages <- ceiling(length(plot_list) / plots_per_page)

# Open PDF device
pdf(OUTPUT_FIGURE, width = 10, height = 10)

# Loop through each page
for (page in 1:num_pages) {
        # Determine the index of the plots for the current page
        plot_indices <- 
                ((page - 1) * plots_per_page + 1):min(page * plots_per_page, length(plot_list))
        # Extract the plots for the current page
        page_plots <- plot_list[plot_indices]
        # Combine plots for this page
        if (length(page_plots) > 0) {
        combined_plot <- plot_grid(plotlist = page_plots, ncol = 2, nrow = 3)
        print(combined_plot)
        }
        }
dev.off()
