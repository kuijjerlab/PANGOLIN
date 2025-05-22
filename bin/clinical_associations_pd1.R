
#####################
## Load R Packages ##
#####################
# List of required libraries
required_libraries <- c("data.table", "dplyr", "doParallel", 
                        "survminer", "survival", "magrittr", 
                        "gridExtra", "ggplot2", "ggrepel",
                        "smplot2", "effectsize", "rstatix",
                        "optparse")

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
            type = "character", default = NULL,
            help = "Path to the cox results for the PD1 pathway",
            metavar = "character"),
       make_option(c("-p", "--tumor_main_dir"), type = "character", 
            default = NULL,
            help = "Path to the tumor main directory",
            metavar = "character"),
        make_option(c("-r", "--results_pd1_groups"), type = "character", 
            default = NULL,
            help = "Path to the result file comparing clinical groups",
            metavar = "character"),
       make_option(c("-s", "--results_pd1_numeric"), type = "character", 
            default = NULL,
            help = "Path to the result file with numerical features",
            metavar = "character")
)

# Parse the command-line arguments
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Assign parsed arguments to variables
COX_RESULTS_FILE <- opt$cox_results_file
TUMOR_MAIN_DIR <- opt$tumor_main_dir
RESULTS_PD1_GROUPS <- opt$results_pd1_groups
RESULTS_PD1_NUMERIC <- opt$results_pd1_numeric

########################
## Load Helper Scripts ##
########################
# source required functions

source("bin/clinical_association_pd1_fn.R")
source("bin/cox_regression_tumor_fn.R")
source("bin/extract_clinical_data_fn.R")
set.seed(1234)
df <- fread(COX_RESULTS_FILE)
pfi_cancer <-  c("brca", "lgg", "prad", "read", "tgct", "thca", "thym")
pfi_cancer <- toupper(pfi_cancer)
# Filter the results based on the type and cancer
df <- df %>%
            filter(type == "PFI" & cancer %in% pfi_cancer |
                    type == "OS" & !cancer %in% pfi_cancer)

df_clean <- df[df$pvalue < 0.05,]
dirs <- list.dirs(TUMOR_MAIN_DIR, recursive = TRUE, full.names = TRUE)
pd1_dirs <- dirs[grepl("pd1_data", basename(dirs))]
clin_dirs <- dirs[grepl("clinical", basename(dirs))]
res_groups_all <- NULL
res_numeric_all <- NULL
for (i in 1:length(df_clean)){
    cancer <- df_clean$cancer[i]
    component <- df_clean$component[i]
    tumor_pd1_dir <- pd1_dirs[grep(cancer, pd1_dirs)]
    tumor_clin_dir <- clin_dirs[grep(cancer, clin_dirs)]
    tumor_clin_file <- 
        list.files(tumor_clin_dir, recursive = TRUE, full.names = TRUE)
    res_groups <- clin_association_groups_pd1(tumor = cancer, 
                        clin_cancer_file =  tumor_clin_file,
                        pd1_dir = tumor_pd1_dir,
                        component = component)
    res_groups$cancer <- cancer
    res_numeric <- clin_association_numeric_pd1(tumor = cancer, 
                        clin_cancer_file =  tumor_clin_file,
                        pd1_dir = tumor_pd1_dir,
                        component = component,
                        correlation_type = c("spearman"))
    res_groups_all <- rbind(res_groups, res_groups_all)  
    res_numeric_all <- rbind(res_numeric, res_numeric_all) 
}
write.table(res_groups_all,
        RESULTS_PD1_GROUPS,
        col.names = T,
        row.names = F, 
        sep = "\t",
        quote = F)


write.table(res_numeric_all,
        RESULTS_PD1_NUMERIC,
        col.names = T,
        row.names = F, 
        sep = "\t",
        quote = F)
