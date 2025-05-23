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
        c("-c", "--cancer_color_file"),
        type = "character",
        default = NULL,
        help = "Path to the cancer color file",
        metavar = "character")    
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

## Initialize variable
DATASETS_TO_PLOT_COLA_CLUSTERS_FILE <- opt$datasets_to_plot_cola_clusters
CANCER_COLOR_FILE <- opt$cancer_color_file


########################
## Load Helper Scripts ##
########################
# source required functions
source("bin/create_cancer_legend_fn.R")
source("bin/sanky_plots_fn.R")

DATASETS_TO_PLOT_COLA_CLUSTERS_FILE <- "/storage/kuijjerarea/tatiana/PANGOLIN/data_all/cola_consensus_clustering/datasets_to_plot_cola_clusters.RData"
CANCER_COLOR_FILE <- "/storage/kuijjerarea/tatiana/PANGOLIN/data_all/samples/cancer_colors.txt"
UMAP_EXP_FILE <-
UMAP_IND_FILE <- 


load(DATASETS_TO_PLOT_COLA_CLUSTERS_FILE, data <- new.env())
datasets_to_plot <- data$datasets_to_plot
uvm_data <- datasets_to_plot$UVM
uvm_data
prad_data <- datasets_to_plot$PRAD
prad_data
p1 <- sanky_plot_viridis(prad_data, "PRAD")
p2 <- sanky_plot_viridis(uvm_data, "UVM")
ggarrange(p1, p2, ncol = 1)







umap_ind_file <- file.path(int_res_dir, "tsne_ind_all_cancers_primary.txt")
umap_exp_file <- file.path(int_res_dir, "tsne_exp_all_cancers_primary.txt")

p3 <- plot_TSNE_a(umap_exp_file, cancer_color_file = cancer_color_file)
p4 <- plot_TSNE_a(umap_ind_file, cancer_color_file = cancer_color_file)




pdf(file.path(fig_dir, "/test_figure_1.pdf"), width = 8, height = 8)
ggarrange(p3, p1, p4, p2, ncol = 2, nrow = 2, labels = c("A.", "C.", "B.", "D."), 
          font.label = list(size = 16, color = "black", face = "bold"))
dev.off()
pdf(file.path(fig_dir, "/test_figure_1.pdf"), width = 24, height = 8)
ggarrange(p3, p4, p1, p2, ncol = 4, nrow = 1, labels = c("A.", "C.", "B.", "D."), 
          font.label = list(size = 16, color = "black", face = "bold"))
dev.off()