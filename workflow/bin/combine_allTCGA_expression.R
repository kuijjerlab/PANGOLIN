
#####################
## Load R packages ##
#####################
# List of packages to be loaded
required_libraries <- c("TCGAbiolinks", "SummarizedExperiment", 
            "dplyr", "DT", "optparse", "data.table")

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
        c("-d", "--expression_dir"),
        type = "character",
        default = NULL,
        help = "Path to the directory with expression files.",
        metavar = "character"),
    optparse::make_option(
        c("-e", "--combined_expression_file"),
            type = "character",
            default = NULL,
            help = "Path to the output combined expression file.",
            metavar = "character"),
    optparse::make_option(
        c("-g", "--group_file"),
        type = "character",
        default = NULL,
        help = "Path to the output group file.",
        metavar = "character"),
    optparse::make_option(
        c("-f", "--feature_file"),
        type = "character",
        default = NULL,
        help = "Path to the output feature file.",
        metavar = "character")
)

### Parse options ###
opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

### Initialize variable ###

EXPRESSION_DIR <- opt$expression_dir
OUTPUT_EXPRESSION_FILE <- opt$combined_expression_file
GROUP_FILE <- opt$group_file
FEATURE_FILE <- opt$feature_file


### load functions ###
source("workflow/bin/download_gdc_fn.R")
exp_files <- list.files(EXPRESSION_DIR, pattern = ".RData")
exp_files

### combine summarized experiments in one matrix ###
exp <- NULL
exp_all <- NULL
groups <- NULL
groups_all <- NULL

for (i in 1:length(exp_files)) {
        project_id <- gsub(".RData", "", exp_files[i])
        exp <- readSummarizedExperiment(file.path(EXPRESSION_DIR, exp_files[i]))
        exp_all <- cbind(exp, exp_all)
        samples <- colnames(exp)
        groups <- data.table("sample_id" = samples, "project" = project_id)
        groups_all <- rbind(groups, groups_all)
}
write.table(exp_all,
            file = OUTPUT_EXPRESSION_FILE,
            col.names = T, row.names = T, sep = "\t", quote = F)


write.table(groups_all, GROUP_FILE,
            col.names = F, row.names = F, sep = "\t", quote = FALSE)

load(file.path(EXPRESSION_DIR, exp_files[1]), data <- new.env())
data <- data[["mrna_df"]]
features <- rowRanges(data)
save(features, file = FEATURE_FILE)