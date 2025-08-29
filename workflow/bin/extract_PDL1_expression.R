#####################
## Load R packages ##
#####################

required_libraries <- c("data.table", "optparse", "tidyr", "dplyr")
for (lib in required_libraries) {
  suppressPackageStartupMessages(library(lib, character.only = TRUE,
                                quietly = TRUE))
}
####################
## Read arguments ##
####################
option_list = list(
    make_option(c("-e", "--exp_file"), type = "character", default = NULL,
            help = "Path to the expression file",
            metavar = "character"),
    make_option(c("-s", "--samples_file"), type = "character", default = NULL,
            help = "Path to the samples file",
            metavar = "character"),
    make_option(c("-g", "--gene_id"), type = "character", default = NULL,
            help = "Gene id to filter (default: CD274)",
            metavar = "character"),
    make_option(c("-t", "--tumor"), type = "character", default = NULL,
            help = "Tumor type to filter (e.g., ACC, BRCA, LUNG)",
            metavar = "character"),
    make_option(c("-o", "--output"), type = "character", default = NULL,
            help = "Output file path",
            metavar = "character")
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)



EXPRESSION_FILE <- opt$exp_file
SAMPLES_FILE <- opt$samples_file
GENE_ID <- opt$gene_id
TUMOR_TYPE <- opt$tumor
OUTPUT_FILE <- opt$output

######################
## Load functions ##
######################
source("workflow/bin/cox_regression_tumor_fn.R")
source("workflow/bin/extract_PDL1_expression_fn.R")

##############################
## Filter for tumor type ##
##############################
exp_pdl1 <- process_gene_expression(exp_file = EXPRESSION_FILE,
                                samples_file = SAMPLES_FILE,
                                gene_id = GENE_ID,
                                tumor = TUMOR_TYPE)

##########################
## Save filtered data ##
##########################
dir.create(dirname(OUTPUT_FILE), recursive = TRUE, showWarnings = FALSE)
fwrite(exp_pdl1, file = OUTPUT_FILE, sep = "\t",
    row.names = FALSE, quote = FALSE, col.names = TRUE)
