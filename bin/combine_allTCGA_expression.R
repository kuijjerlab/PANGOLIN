
#####################
## Load R packages ##
#####################
# List of packages to be loaded
required_libraries <- c("TCGAbiolinks", "SummarizedExperiment", 
            "dplyr", "DT", "optparse")

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
        c("-d", "--directory_expression"),
        type = "character",
        default = NULL,
        help = "Path to the directory with expression files.",
        metavar = "character"),
    optparse::make_option(
        c("-o", "--output_dir"),
        type = "character",
        default = NULL,
        help = "Path to the the output directory.",
        metavar = "character"))

### Parse options ###
opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

### Initialize variable ###

EXPRESSION_DIR <- opt$directory_expression
OUTPUT_DIR <- opt$output_dir

### load functions ###
source("bin/download_gdc_fn.R")
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
save(exp_all, file = file.path(OUTPUT_DIR, "hg38_STAR_counts.RData"))

write.table(groups_all, file.path(OUTPUT_DIR, "hg38_sample_groups.tsv"),
            col.names = F, row.names = F, sep = "\t", quote = FALSE)


load(file.path(exp_files_dir, exp_files[1]), data <- new.env())
data <- data[["mrna_df"]]
features <- rowRanges(data)
save(features, file = file.path(OUTPUT_DIR, "hg38_features.RData"))
