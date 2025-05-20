library(ComplexHeatmap)
library(data.table)
library(ggplot2)
library(circlize)
library(cowplot)
library(grid)
library(plyr)
library(dplyr)

#####################
## Load R Packages ##
#####################
# List of required libraries
required_libraries <- c("data.table", "optparse", "ComplexHeatmap", "ggplot2",
                        "circlize", "cowplot", "grid", "plyr", "dplyr")
# Load each library and suppress startup messages
for (lib in required_libraries) {
  suppressPackageStartupMessages(library(lib, character.only = TRUE,
                                quietly = TRUE))
}

###################
## Parse Arguments ##
####################
# Define command-line options
option_list <- list(
    optparse::make_option(
        c("-p", "--pcp_results_all_cancers_file"),
        type = "character",
        default = NULL,
        help = "Path to the porcupine results file (all cancers).",
        metavar = "character"),
optparse::make_option(
        c("-f", "--pathways_hierarchy_file"),
        type = "character",
        default = NULL,
        help = "Path to the pathways hierarchy file.",
        metavar = "character"),
optparse::make_option(
        c("-m", "--pathways_hsa_id_file"),
        type = "character",
        default = NULL,
        help = "Path to the pathways_hsa_id_file (reactome name - hsa id).",
        metavar = "character"),
optparse::make_option(
        c("-l", "--list_of_pathways_file"),
        type = "character",
        default = NULL,
        help = "Path to the list of pathways file.",
        metavar = "character"), 
optparse::make_option(
        c("-g", "--figure_pathway_intersection"),
        type = "character",
        default = NULL,
        help = "Path to the output figure file with number of shared porcupine
        pathways between cancers.",
        metavar = "character"),
optparse::make_option(
        c("-f", "--figure_shared_categories"),
        type = "character",
        default = NULL,
        help = "Path to the output figure file with shared categories.",
        metavar = "character")
)

# Parse the command-line arguments
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Assign parsed arguments to variables
## Initialize variable
PCP_RESULTS_FILE = opt$pcp_results_all_cancers_file
PATHWAYS_HIERARCHY_FILE = opt$pathways_hierachy_file
PATHWAYS_HSA_ID_FILE = opt$pathways_hsa_id_file
LIST_PATHWAYS_FILE = opt$list_of_pathways_file
FIG_PATHWAY_INTERSECTION = opt$figure_pathway_intersection
FIG_SHARED_CATEGORIES = opt$figure_shared_categories

########################
## Load Helper Scripts ##
########################
# source required functions
set.seed(1234)
source("bin/analyze_pcp_results_fn.R")
source("bin/plot_PORCUPINE_results_fn.R")

# Assign biological groups to the filtered PORCUPINE pathways
res_all <- fread(PCP_RESULTS_FILE)
head(res_all)
ptws_sel <- data.table("pathway" = unique(res_all$pathway))
ptws_dat <- assign_functions_to_pathways(
                        pathways_hsa_ids_file = PATHWAYS_HSA_ID_FILE,
                        pathways_hierachy_file = PATHWAYS_HIERARCHY_FILE,
                        list_of_pathways_file = LIST_PATHWAYS_FILE, 
                        pathways = ptws_sel)


n_pathways <- table(res_all$cancer)
res_list <- split(res_all$pathway, res_all$cancer)
intersection_matrix <- intersect_pathways(res_list, type = c("jaccard"))
colnames(intersection_matrix) <- toupper(colnames(intersection_matrix))
rownames(intersection_matrix) <- toupper(rownames(intersection_matrix))

# plot the intersection results
pdf(FIG_PATHWAY_INTERSECTION, width = 10, height = 10)
plot_intersection_heatmap(
            intersection_matrix = intersection_matrix, 
            pathways_counts = n_pathways)
dev.off()

res_pathways <- 
    preprocess_pathway_results(res_all = res_all, ptws_dat = ptws_dat)

categories_shared <- res_pathways$categories_shared
ptws_shared <- res_pathways$ptws_shared
stat <- res_pathways$stat
col_categories <- res_pathways$col_categories
dotplot <- plot_pathway_category_dotplot(stat = stat,
                                    col_categories = col_categories)
ht_list <- make_category_heatmap_list(categories_shared = categories_shared,
                                    ptws_shared = ptws_shared,
                                    col_categories =  col_categories,
                                    intersection_matrix = intersection_matrix)
ht_list_up <- Reduce(`%v%`, ht_list)
final_plot <- combine_dotplot_and_heatmap(dotplot, ht_list_up)
pdf(FIG_SHARED_CATEGORIES, width = 12, height =12)
print(final_plot)
dev.off()
