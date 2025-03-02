#####################
## Load R packages ##
#####################
required_libraries <- c("data.table", "dplyr", "doParallel",
            "survminer", "survival", "magrittr", "gridExtra",
            "purrr", "optparse")
for (lib in required_libraries) {
  suppressPackageStartupMessages(library(lib, character.only = TRUE,
                                quietly = TRUE))
}

####################
## Read arguments ##
####################
### Options
options(stringsAsFactors = FALSE)

### Command line options
option_list <- list(
    optparse::make_option(
        c("-t", "--tumor_clin_file_path"),
        type = "character",
        default = NULL,
        help = "Path to the clinical file.",
        metavar = "character"),
    optparse::make_option(
        c("-p", "--tumor_pd1_dir"),
        type = "character",
        default = NULL,
        help = "Path to the pd1 tumor directory.",
        metavar = "character"),
    optparse::make_option(
        c("-f", "--number_folds"),
        type = "integer",
        default = 5,
        help = "Number of folds for cross-validation.",
        metavar = "integer"),
    optparse::make_option(
        c("-m", "--number_times"),
        type = "integer",
        default = 10,
        help = "Number of times to run a model.",
        metavar = "integer"),    
    optparse::make_option(
        c("-c", "--number_cores"),
        type = "integer",
        default = 1,
        help = "Number of cores.",
        metavar = "integer"),
    optparse::make_option(
        c("-a", "--alpha"),
        type = "numeric",
        default = 1,
        help = "Alpha value for the elastic net.",
        metavar = "numeric"),
    optparse::make_option(
        c("-o", "--output_file"),
        type = "character",
        default = NULL,
        help = "Path to the output_file",
        metavar = "character"))

opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

## Initialize variable
TUMOR_CLIN_FILE <- opt$tumor_clin_file_path
TUMOR_PD1_DIR <- opt$tumor_pd1_dir
NUMBER_FOLDS <- opt$number_folds
NUMBER_TIMES <- opt$number_times
NUMBER_CORES <- opt$number_cores
ALPHA <- opt$alpha
OUTPUT_FILE <- opt$output_file

source("bin/cox_regression_tumor_fn.R")
dir.create(dirname(OUTPUT_FILE), recursive = TRUE, showWarnings = FALSE)
res <- run_glmnet_ntimes_pd1(tumor_clin_file_path = TUMOR_CLIN_FILE,
                    tumor_pd1_dir = TUMOR_PD1_DIR,
                    number_folds = NUMBER_FOLDS,
                    ntimes = NUMBER_TIMES,
                    ncores = NUMBER_CORES, 
                    alpha = ALPHA)
write.table(res, OUTPUT_FILE, col.names = TRUE,
    sep = "\t", row.names = FALSE, quote = FALSE)



TUMOR_CLIN_FILE <- "data_individual_cancers/BRCA/clinical/curated_clinical_BRCA.txt"
TUMOR_PD1_DIR <- "data_individual_cancers/BRCA/pd1_data"
NUMBER_FOLDS <- 10
NUMBER_CORES <- 5
ALPHA <- 0.3
NUMBER_TIMES <- 10
OUTPUT_FILE <- "data_individual_cancers/BRCA/cox/BRCA_cox_multivariate_res.txt"



data <- combine_info_for_cancer(TUMOR_CLIN_FILE, 
                               TUMOR_PD1_DIR)

        
pd1_links <- load_pd1_generic( TUMOR_PD1_DIR, type = "pd1_links")
pd1_net <- load_pd1_generic( TUMOR_PD1_DIR, type = "pd1_net")
            rownames(pd1_net) <- pd1_links
            data$pd1_net <- pd1_net
            data$pd1_scores <- load_pd1_generic(tumor_pd1_dir,
                                type = "pd1_scores")


data <- load_pd1_object("/storage/kuijjerarea/tatiana/PANGOLIN/data_individual_cancers/BRCA/pd1_data/pd1_net_norm_BRCA.RData", object_name = "pd1_net")
data