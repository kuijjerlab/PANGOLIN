#####################
## Load R packages ##
#####################
required_libraries <- c("TCGAbiolinks", "SummarizedExperiment", 
                        "dplyr", "DT", "optparse", "data.table")

for (lib in required_libraries) {
    suppressPackageStartupMessages(library(lib, character.only = TRUE,
                                  quietly = TRUE))
}

####################
## Read arguments ##
####################
options(stringsAsFactors = FALSE)

option_list <- list(
    optparse::make_option(
        c("-f", "--gdc_files"),
        type = "character",
        default = NULL,
        help = "Comma or space-separated list of GDC RData files.",
        metavar = "character"
    ),
    optparse::make_option(
        c("-e", "--combined_expression_file"),
        type = "character",
        default = NULL,
        help = "Path to the output combined expression file.",
        metavar = "character"
    ),
    optparse::make_option(
        c("-g", "--group_file"),
        type = "character",
        default = NULL,
        help = "Path to the output group file.",
        metavar = "character"
    ),
    optparse::make_option(
        c("-r", "--feature_file"),
        type = "character",
        default = NULL,
        help = "Path to the output feature file.",
        metavar = "character"
    )
)

opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

# Parse file list
GDC_FILES <- opt$gdc_files
OUTPUT_EXPRESSION_FILE <- opt$combined_expression_file
GROUP_FILE <- opt$group_file
FEATURE_FILE <- opt$feature_file

# Convert space-separated string to vector
exp_files <- unlist(strsplit(GDC_FILES, " "))
exp_files <- exp_files[exp_files != ""]  # Remove empty strings

### load functions ###
source("workflow/bin/download_gdc_fn.R")

### combine summarized experiments in one matrix ###
exp_all <- NULL
groups_all <- NULL

for (i in 1:length(exp_files)) {
    project_id <- gsub(".RData", "", basename(exp_files[i]))
    exp <- readSummarizedExperiment(exp_files[i])  # Use full path
    exp_all <- cbind(exp, exp_all)
    samples <- colnames(exp)
    groups <- data.table("sample_id" = samples, "project" = project_id)
    groups_all <- rbind(groups, groups_all)
}

write.table(exp_all,
            file = OUTPUT_EXPRESSION_FILE,
            col.names = TRUE, row.names = TRUE, sep = "\t", quote = FALSE)

write.table(groups_all, GROUP_FILE,
            col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)

load(exp_files[1], data <- new.env())
data <- data[["mrna_df"]]
features <- rowRanges(data)
save(features, file = FEATURE_FILE)