#####################
## Load R packages ##
#####################
required_libraries <- c("data.table", "ggplot2", "ggforce")
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
        c("-b", "--batch_results_dir"),
        type = "character",
        default = NULL,
        help = "Path to the batch results directory.",
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

BATCH_DIR <- opt$batch_results_dir
OUTPUT_FILE <- opt$output_file

# source the function to combine results
source("bin/analyze_batch_fn.R")

# combine and plot the MBatch results
cancers <- list.files(BATCH_DIR)
cancers <- rev(cancers)

res_all <- NULL
res <- NULL
for (m in 1:length(cancers)){
  cancer <- cancers[m]
  res <- combine_mbatch_results(BATCH_DIR, cancer)
  res <- res[grep("1", res$Annotation), ]
  res_all <- rbind(res_all, res)
}

plot_mbatch_dsc(res_all, output_file = OUTPUT_FILE)
