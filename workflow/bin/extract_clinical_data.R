#####################
## Load R packages ##
#####################
required_libraries <- c("data.table", "openxlsx", "optparse")
for (lib in required_libraries) {
  suppressPackageStartupMessages(library(lib, character.only = TRUE,
                                quietly = TRUE))
}


####################
## Read arguments ##
####################
option_list = list(
    make_option(c("-c", "--clinical"), type = "character", default = NULL,
            help = "Path to the clinical data Excel file",
            metavar = "character"),
    make_option(c("-t", "--tumor"), type = "character", default = NULL,
            help = "Tumor type to filter (e.g., ACC, BRCA, LUNG)",
            metavar = "character"),
    make_option(c("-o", "--output"), type = "character", default = NULL,
            help = "Output file path for the filtered clinical data",
            metavar = "character")
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)



CLINICAL_FILE <- opt$clinical
TUMOR_TYPE <- opt$tumor
OUTPUT_FILE <- opt$output

message("Processing clinical data for tumor type: ", TUMOR_TYPE)

######################
## Load functions ##
######################
source("workflow/bin/extract_clinical_data_fn.R")

######################
## Load clinical data ##
######################
clinical_data <- load_clin_curated(CLINICAL_FILE)

##############################
## Filter for tumor type ##
##############################
filtered_data <- select_tumor_clin_curated(clinical_data, TUMOR_TYPE)

##########################
## Save filtered data ##
##########################
dir.create(dirname(OUTPUT_FILE), recursive = TRUE, showWarnings = FALSE)
fwrite(filtered_data, file = OUTPUT_FILE, sep = "\t",
    row.names = FALSE, quote = FALSE)


