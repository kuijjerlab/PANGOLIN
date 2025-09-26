
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
        make_option(c("-g", "--combined_results_groups"),
                type = "character",
                default = NULL,
                help = "Path to the combined categorical results file",
                metavar = "character"),
        make_option(c("-n", "--combined_results_numeric"),
                type = "character",
                default = NULL,
                help = "Path to the combined numeric results file",
                metavar = "character"),
        make_option(c("-p", "--pd1_scores_files"),
                type = "character", 
                default = NULL,
                help = "Comma-separated list of PD1 scores files",
                metavar = "character"),
        make_option(c("-l", "--clinical_files"),
                type = "character",
                default = NULL,
                help = "Comma-separated list of clinical files",
                metavar = "character"),
        make_option(c("-t", "--cancer_types"),
                type = "character",
                default = NULL,
                help = "Comma-separated list of cancer types",
                metavar = "character"),
        make_option(c("-o", "--output_figure_file"),
                type = "character",
                default = NULL,
                help = "Path to the output figure file",
                metavar = "character"),
        make_option(c("-a", "--padjust_threshold"),
                type = "numeric",
                default = 0.01,
                help = "P-adjusted threshold for significance (default: 0.01)",
                metavar = "numeric"),
        make_option(c("-e", "--effect_size_threshold"),
                type = "numeric",
                default = 0.4,
                help = "Effect size threshold (default: 0.4)",
                metavar = "numeric")
)

# Parse the command-line arguments
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Assign parsed arguments to variables
COX_RESULTS_FILE <- opt$cox_results_file
RESULTS_CATEGORICAL_FILE <- opt$combined_results_groups
RESULTS_NUMERIC_FILE <- opt$combined_results_numeric
PD1_SCORES_FILES <- unlist(strsplit(opt$pd1_scores_files, ","))
CLINICAL_FILES <- unlist(strsplit(opt$clinical_files, ","))
CANCER_TYPES <- unlist(strsplit(opt$cancer_types, ","))
OUTPUT_FIGURE <- opt$output_figure_file
PADJUST_THRESHOLD <- opt$padjust_threshold
EFFECT_SIZE_THRESHOLD <- opt$effect_size_threshold

# Validate inputs
if (length(PD1_SCORES_FILES) != length(CLINICAL_FILES) ||
    length(PD1_SCORES_FILES) != length(CANCER_TYPES)) {
    stop("Number of PD1 scores files, clinical files, and cancer types must match")
}

cat("Using p-adjusted threshold:", PADJUST_THRESHOLD, "\n")
cat("Using effect size threshold:", EFFECT_SIZE_THRESHOLD, "\n")



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
      filter(padjust <= PADJUST_THRESHOLD) %>%
      filter(abs(effect_size_cor) >= EFFECT_SIZE_THRESHOLD)

res_pc1 <- data_all %>%
      filter(principal_component == "PC1") %>%
      mutate(feature_cancer = paste0(cancer, "_", clin_feature))


features_pc1 <- unique(res_pc1$feature_cancer)
res_pc2 <- data_all %>%
      filter(principal_component == "PC2") %>%
      mutate(feature_cancer = paste0(cancer, "_", clin_feature))


features_pc2 <- unique(res_pc2$feature_cancer)

# Create lookup tables for files by cancer type
pd1_lookup <- setNames(PD1_SCORES_FILES, toupper(CANCER_TYPES))
clinical_lookup <- setNames(CLINICAL_FILES, toupper(CANCER_TYPES))

plots_pc1 <- list()

for (i in seq_along(features_pc1)) {
      splits <- split_string_function(features_pc1[i])
      tumor <- toupper(splits$cancer)
      feature_to_plot <- splits$feature
      
      # Get files for this cancer type
      tumor_clin_file <- clinical_lookup[[tumor]]
      pd1_scores_file <- pd1_lookup[[tumor]]
      
      if (is.na(tumor_clin_file) || is.na(pd1_scores_file)) {
            cat("Warning: Missing files for cancer type", tumor, "\n")
            next
      }
      
      plot <- plot_clin_feature(tumor = tumor,
                              results_categorical = results_categorical_all,
                              feature_to_plot = feature_to_plot,
                              component = "PC1",
                              clin_cancer_file = tumor_clin_file,
                              pd1_scores_file = pd1_scores_file)
      plots_pc1[[i]] <- plot
}

plots_pc2 <- list()

for (i in seq_along(features_pc2)) {
      splits <- split_string_function(features_pc2[i])
      tumor <- toupper(splits$cancer)
      feature_to_plot <- splits$feature
      
      # Get files for this cancer type
      tumor_clin_file <- clinical_lookup[[tumor]]
      pd1_scores_file <- pd1_lookup[[tumor]]
      
      if (is.na(tumor_clin_file) || is.na(pd1_scores_file)) {
            cat("Warning: Missing files for cancer type", tumor, "\n")
            next
      }
      
      plot <- plot_clin_feature(tumor = tumor,
                              results_categorical = results_categorical_all,
                              feature_to_plot = feature_to_plot,
                              component = "PC2",
                              clin_cancer_file = tumor_clin_file,
                              pd1_scores_file = pd1_scores_file)
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
