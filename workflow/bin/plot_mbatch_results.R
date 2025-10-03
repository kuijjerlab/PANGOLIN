#####################
## Load R packages ##
#####################
required_libraries <- c("data.table", "ggplot2", "ggforce", "optparse")
for (lib in required_libraries) {
    suppressPackageStartupMessages(
        library(lib, character.only = TRUE, quietly = TRUE)
    )
}
####################
## Read arguments ##
####################
option_list <- list(
    optparse::make_option(
        c("-b", "--batch_files"),
        type = "character",
        default = NULL,
        help = "Path to the batch files.",
        metavar = "character"
    ),
    optparse::make_option(
        c("-o", "--output_file"),
        type = "character",
        default = NULL,
        help = "Path to the output file to save the plot.",
        metavar = "character"
    )
)
opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

### Initialize variable ###

BATCH_FILES <-  unlist(strsplit(opt$batch_files, " "))
OUTPUT_FILE <- opt$output_file

# source the function to combine results
source("workflow/bin/analyze_batch_fn.R")

# combine and plot the MBatch results

res <- combine_mbatch_results(BATCH_FILES)
res <- res[-grep("Purity_singscore", res$batch), ]  # remove purity singscore results
# Filter for DSC (1,2) annotations specifically
res_all <- res[grep("1", res$Annotation), ]
plot_mbatch_dsc(res_all, output_file = OUTPUT_FILE)
