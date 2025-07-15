# Use the environment variable set by Snakemake
# Auto-detect conda environment (Snakemake will set this)
python_env <- Sys.getenv("CONDA_PREFIX")
if (python_env != "" && python_env != "None") {
    Sys.setenv(MBATCH_PYTHON_ENV = python_env)
    cat("Using Snakemake conda environment:", python_env, "\n")
} else {
    stop("No conda environment detected. Make sure Snakemake is using --use-conda")
}

# Sys.setenv(MBATCH_PYTHON_ENV = "/storage/kuijjerarea/tatiana/anaconda3/envs/mbatch_minimal")

# require package installation
# devtools::install_github("MD-Anderson-Bioinformatics/BatchEffectsPackage/apps/MBatch", force = TRUE)

#####################
## Load R packages ##
#####################
# List of packages to be loaded
required_libraries <- c("MBatch", "data.table", "optparse")

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
        c("-t", "--tumor_type"),
        type = "character",
        default = NULL,
        help = "Tumor type.",
        metavar = "character"),
    optparse::make_option(
        c("-e", "--expression_file"),
        type = "character",
        default = NULL,
        help = "Path to the expression file (cancer specific).",
        metavar = "character"),
    optparse::make_option(
        c("-g", "--group_file"),
        type = "character",
        default = NULL,
        help = "Path to the group file.",
        metavar = "character"),
    optparse::make_option(
        c("-f", "--feature_file"),
        type = "character",
        default = NULL,
        help = "Path to the feature file.",
        metavar = "character"),
    optparse::make_option(
        c("-b", "--batch_file"),
        type = "character",
        default = NULL,
        help = "Path to the batch file.",
        metavar = "character"),
    optparse::make_option(
        c("-c", "--clin_file"),
        type = "character",
        default = NULL,
        help = "Path to the clinical file.",
        metavar = "character"),
    optparse::make_option(
        c("-o", "--output_directory"),
        type = "character",
        default = NULL,
        help = "Path to the output directory (cancer specific).")
)

### Parse options ###
opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

### Initialize variable ###

TUMOR_TYPE <- opt$tumor_type
EXPRESSION_FILE <- opt$expression_file
GROUP_FILE <- opt$group_file
FEATURE_FILE <- opt$feature_file
BATCH_FILE <- opt$batch_file
CLIN_FILE <- opt$clin_file
OUTPUT_DIR <- opt$output_directory


### source functions ###
source("bin/analyze_batch_fn.R")

cat("Loading expression file...\n")
load(EXPRESSION_FILE, exp <- new.env())
log2exp <- exp[['exp_proj']]

cat("Loading other data files...\n")
load(FEATURE_FILE, features <- new.env())
features <- features[['features']]

groups <- fread(GROUP_FILE, head = F)
batch_info <- load_batch(BATCH_FILE)
batch_info <- batch_info[batch_info$Tissues %in% c("cancer")]
batch_info$Tissues <- NULL
batch_info[is.na(batch_info)] <- c("not_available")
clin <- load_clin_rdata(CLIN_FILE)
batch_info$platform <- 
    clin$gdc_platform[match(batch_info$Samples, 
    clin$gdc_cases.samples.portions.analytes.aliquots.submitter_id)]
batch_info[is.na(batch_info)] <- c("not_available")

# Create the output directory if it doesn't exist
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Set working directory to the output directory
setwd(OUTPUT_DIR)

batch_process_cancer(tumor_type = TUMOR_TYPE,
                                 log2exp = log2exp,
                                 groups = groups, 
                                 batch_info = batch_info, 
                                 clin = clin)



