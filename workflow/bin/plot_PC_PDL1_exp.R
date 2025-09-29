#####################
## Load R packages ##
#####################
required_libraries <- c("data.table", "dplyr", "stringr",
            "ggplot2", "cowplot", "optparse", "smplot2")
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
        c("-c", "--cox_summary_all_cancers"),
        type = "character",
        default = NULL,
        help = "Path to the cox summary file.",
        metavar = "character"),
    optparse::make_option(
        c("-d", "--combined_patient_data_files"),
        type = "character",
        default = NULL,
        help = "Path to the combined patient data files.",
        metavar = "character"),
    optparse::make_option(
        c("-o", "--output_file"),
        type = "character",
        default = NULL,
        help = "Path to the output pdf figure file.",
        metavar = "character")
        )
opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

## Initialize variable
COX_SUMMARY_ALL <- opt$cox_summary_all_cancers
COMBINED_PATIENT_DATA_FILES <- opt$combined_patient_data_files
OUTPUT_PDF_FILE <- opt$output_file

source("workflow/bin/merge_patient_data_fn.R")
source("workflow/bin/plotting_fn.R")

dir.create(dirname(OUTPUT_PDF_FILE), recursive = TRUE, showWarnings = FALSE)


COMBINED_PATIENT_DATA_FILES <- 
    unlist(strsplit(opt$combined_patient_data_files, " "))
plot_list <- generate_pc_cd274_plots(cox_results_file = COX_SUMMARY_ALL,
                            combined_patient_data_files = COMBINED_PATIENT_DATA_FILES) 

# Determine the number of plots per page and calculate number of pages needed
plots_per_page <- 16
num_pages <- ceiling(length(plot_list) / plots_per_page)

# Open PDF device
pdf(OUTPUT_PDF_FILE, width = 12, height = 12)

# Loop through each page
for (page in 1:num_pages) {
        # Determine the index of the plots for the current page
        plot_indices <- 
                ((page - 1) * plots_per_page + 1):min(page * plots_per_page, length(plot_list))
        # Extract the plots for the current page
        page_plots <- plot_list[plot_indices]
        # Combine plots for this page
        if (length(page_plots) > 0) {
        combined_plot <- plot_grid(plotlist = page_plots, ncol = 4)
        print(combined_plot)
        }
    }
dev.off()
