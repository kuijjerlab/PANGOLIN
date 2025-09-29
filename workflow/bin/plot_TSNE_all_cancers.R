#####################
## Load R Packages ##
#####################
# List of required libraries
required_libraries <- c("data.table", "ggplot2", "optparse",
                        "ggsankey", "ggpubr")
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
    make_option(
        c("-z", "--datasets_to_plot_cola_clusters"),
        type = "character",
        default = NULL,
        help = "Path to a file containing the information about the 
                datasets to plot for cola cluster comparison",
        metavar = "character"),
    make_option(
        c("-f", "--tsne_data_expression"),
        type = "character",
        default = NULL,
        help = "Path to the file with tsne results for expression.",
        metavar = "character"),
    make_option(
        c("-m", "--tsne_data_indegree"),
        type = "character",
        default = NULL,
        help = "Path to the file with tsne results for indegree.",
        metavar = "character"),
    make_option(
        c("-c", "--cancer_color_file"),
        type = "character",
        default = NULL,
        help = "Path to the cancer color file",
        metavar = "character"),
    make_option(
        c("-o", "--output_figure_file"),
        type = "character",
        default = NULL,
        help = "Path to the output figure file",
        metavar = "character")
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

## Initialize variable
DATASETS_TO_PLOT_COLA_CLUSTERS_FILE <- opt$datasets_to_plot_cola_clusters
CANCER_COLOR_FILE <- opt$cancer_color_file
TNSE_EXPRESSION_FILE <- opt$tsne_data_expression
TSNE_INDEGREE_FILE <-  opt$tsne_data_indegree
OUTPUT_FIGURE_FILE <- opt$output_figure_file

########################
## Load Helper Scripts ##
########################
# source required functions
source("workflow/bin/create_cancer_legend_fn.R")
source("workflow/bin/sanky_plots_fn.R")

# loading the datasets to plot for cola cluster comparison
load(DATASETS_TO_PLOT_COLA_CLUSTERS_FILE, data <- new.env())
datasets_to_plot <- data$datasets_to_plot
uvm_data <- datasets_to_plot$UVM
uvm_data
prad_data <- datasets_to_plot$PRAD
prad_data
p1 <- sanky_plot_viridis(prad_data, "PRAD")
p2 <- sanky_plot_viridis(uvm_data, "UVM")
ggarrange(p1, p2, ncol = 1)

p3 <- plot_TSNE_all_cancers(tsne_res_file = TNSE_EXPRESSION_FILE,
                            cancer_color_file = CANCER_COLOR_FILE)
p4 <- plot_TSNE_all_cancers(tsne_res_file = TSNE_INDEGREE_FILE, 
                            cancer_color_file = CANCER_COLOR_FILE)

spacer <- ggplot() + theme_void()
pdf(OUTPUT_FIGURE_FILE, width = 10, height = 10)
ggarrange(p3, p1, p4, p2, 
        spacer, spacer,
        ncol = 2, nrow = 3, labels = c("A.", "C.", "B.", "D."), 
        font.label = list(size = 16, color = "black", face = "bold"),
        heights = c(1, 1, 0.3))
dev.off()
