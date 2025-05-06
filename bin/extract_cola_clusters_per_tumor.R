#####################
## Load R packages ##
#####################
required_libraries <- c("data.table", "openxlsx", "optparse")
for (lib in required_libraries) {
  suppressPackageStartupMessages(library(lib, character.only = TRUE,
                                quietly = TRUE))
}


####################
## Read arguments ##
####################
option_list = list(
    make_option(c("-c", "--cluster_file_expression"), type = "character", 
                default = NULL,
                help = "Path to the cluster file expression",
                metavar = "character"),
    make_option(c("-m", "--cluster_file_indegree"), type = "character",
                 default = NULL,
                help = "Path to the cluster file indegree",
                metavar = "character"),
    make_option(c("-t", "--tumor"), type = "character", default = NULL,
                help = "Tumor type to filter (e.g., ACC, BRCA, LUNG)",
                metavar = "character"),
    make_option(c("-e", "--cluster_expression_per_tumor"), type = "character", 
                default = NULL,
                help = "Output file path for the filtered
                clustered data for a specific cancer (expression)",
                metavar = "character"),
    make_option(c("-i", "--cluster_indegree_per_tumor"), type = "character",
                default = NULL,
                help = "Output file path for the filtered 
                clustered data for a specific cancer (indegree)",
                metavar = "character")
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)



CLUSTER_FILE_EXP <- opt$cluster_file_expression
CLUSTER_FILE_IND <- opt$cluster_file_indegree
TUMOR_TYPE <- opt$tumor
EXPRESSION_CLUSTERS_PER_CANCER <- opt$cluster_expression_tumor
INDEGREE_CLUSTERS_PER_CANCER <- opt$cluster_indegree_tumor

source("bin/cola_clustering_fn.R")

##############################
## Filter for tumor type ##
##############################
filtered_data_exp <- 
        filter_cola_clusters_for_tumor(CLUSTER_FILE_EXP, TUMOR_TYPE)
filtered_data_ind <-  
        filter_cola_clusters_for_tumor(CLUSTER_FILE_IND, TUMOR_TYPE)

## Save filtered data ##
##########################
#dir.create(dirname(EXPRESSION_CLUSTERS_PER_CANCER), recursive = TRUE, showWarnings = FALSE)
fwrite(filtered_data_exp, file = EXPRESSION_CLUSTERS_PER_CANCER, sep = "\t",
    row.names = FALSE, quote = FALSE)

#dir.create(dirname(EXPRESSION_CLUSTERS_PER_CANCER), recursive = TRUE, showWarnings = FALSE)
fwrite(filtered_data_ind, file = INDEGREE_CLUSTERS_PER_CANCER, sep = "\t",
    row.names = FALSE, quote = FALSE)