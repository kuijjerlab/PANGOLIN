#####################
## Load R Packages ##
#####################
# List of required libraries
required_libraries <- c("data.table", "dplyr", "optparse",
                        "ggplot2", "ggrepel", "cowplot")
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
        make_option(c("-m", "--cancer_color_file"), 
                type = "character", 
                default = NULL,
                help = "Path to the cancer color file",
                metavar = "character"), 
        make_option(c("-p", "--pc_clinical_association_figure"), 
                type = "character", 
                default = NULL,
                help = "Path to the cancer color file",
                metavar = "character")
)

# Parse the command-line arguments
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Assign parsed arguments to variables
COX_RESULTS_FILE <- opt$cox_results_file
RESULTS_CATEGORICAL_FILE <- opt$results_pd1_groups
RESULTS_NUMERIC_FILE <- opt$results_pd1_numeric
CANCER_COLOR_FILE <- opt$cancer_color_file
PC_CLIN_ASSOCIATIONS_FIGURE <- opt$pc_clinical_association_figure

########################
## Load Helper Scripts ##
########################
# source required functions

source("bin/plot_clinical_associations_fn.R")


res_cat <- process_categorical_results(
                        res_categorical_file = RESULTS_CATEGORICAL_FILE,
                        coxph_results_file = COX_RESULTS_FILE)
res_num <- process_numeric_results(
                        res_numeric_file = RESULTS_NUMERIC_FILE,
                        coxph_results_file = COX_RESULTS_FILE)
coxph_results <- load_coxph_results(coxph_results_file = COX_RESULTS_FILE)
data_all <- combine_data(res_cat, res_num)
colors <- load_cancer_colors(CANCER_COLOR_FILE)
gplot_list <- list()
for (i in 1:nrow(coxph_results)) {
        cancer_PC <- coxph_results$cancer_component[i]
        data <- data_all %>%
                filter(cancer_component == cancer_PC)
        tumor <- unique(data$cancer)
        gplot_list[[i]] <- plot_clinical_associations(data, pc_component,
                        cancer_PC, colors)
}


pdf(PC_CLIN_ASSOCIATIONS_FIGURE, width = 12, height = 12)
combined_plot <- plot_grid(plotlist = gplot_list, ncol = 4)
print(combined_plot)
dev.off()
