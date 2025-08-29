#####################
## Load R packages ##
#####################
required_libraries <- c("data.table", "ggplot2", "optparse", "gridExtra")
for (lib in required_libraries) {
  suppressPackageStartupMessages(library(lib, character.only = TRUE,
                                quietly = TRUE))
}


####################
## Read arguments ##
####################
option_list = list(
    make_option(c("-c", "--cancer_color_file"),
            type = "character", default = NULL,
            help = "Path to the cancer color file",
            metavar = "character"),
    make_option(c("-f", "--figure_dir"), type = "character", default = NULL,
            help = "File path to the output directory",
            metavar = "character")
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

CANCER_COLOR_FILE <- opt$cancer_color_file
FIGURE_DIR <- opt$figure_dir

######################
## Load functions ##
######################
source("workflow/bin/create_cancer_legend_fn.R")


######################
## Load the data ##
######################
legend_data <- fread(CANCER_COLOR_FILE)
legend_data$cancer_type <- toupper(legend_data$cancer_type)
legend_plot <- create_legend(legend_data)
dir.create(dirname(FIGURE_DIR), recursive = TRUE, showWarnings = FALSE)

######################
## Make the plot ##
######################


pdf(file.path(FIGURE_DIR, "cancer_legend.pdf"), width = 10, height = 8)
# Display the legend
gridExtra::grid.arrange(legend_plot)
dev.off()

