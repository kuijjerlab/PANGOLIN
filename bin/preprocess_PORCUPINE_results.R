#####################
## Load R Packages ##
#####################
# List of required libraries
required_libraries <- c("data.table", "optparse")
# Load each library and suppress startup messages
for (lib in required_libraries) {
  suppressPackageStartupMessages(library(lib, character.only = TRUE,
                                quietly = TRUE))
}

# Set a random seed for reproducibility
set.seed(1234)


####################
## Parse Arguments ##
####################
# Define command-line options
option_list <- list(
    optparse::make_option(
        c("-p", "--porcpupine_file_path"),
        type = "character",
        default = NULL,
        help = "Path to the porcupine results file.",
        metavar = "character"),

optparse::make_option(
        c("-f", "--filtered_porcupine_file_path"),
        type = "character",
        default = NULL,
        help = "Path to the filtered porcupine results.",
        metavar = "character")
)

# Parse the command-line arguments
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Assign parsed arguments to variables
## Initialize variable
PCP_FILE <- opt$porcpupine_file_path
PCP_FILT_FILE <- opt$filtered_porcupine_file_path
########################
## Load Helper Scripts ##
########################
# source required functions
source("bin/analyze_pcp_results_fn.R")
res_filt <- filter_pcp_results(PCP_FILE)
write.table(res_filt,
            PCP_FILT_FILE,
            col.names = T,
            row.names = F,
            sep = "\t",
            quote = F)



