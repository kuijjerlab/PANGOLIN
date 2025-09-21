#####################
## Load R packages ##
#####################

required_libraries <- c("data.table", "optparse",
                        "tidyr", "dplyr", "irlba", "Rtsne")
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
    make_option(
        c("-d", "--tumor_dir"),
        type = "character",
        default = NULL,
        help = "Path to the the main tumor directory.",
        metavar = "character"),
    make_option(
        c("-o", "--output_file_expression"),
        type = "character",
        default = NULL,
        help = "Path to the output txt file with tsne results for expression.",
        metavar = "character"),
    make_option(
        c("-m", "--output_file_indegree"),
        type = "character",
        default = NULL,
        help = "Path to the output txt file with tsne results for indegree.",
        metavar = "character")   
        )

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)


## Initialize variable
EXPRESSION_FILE <- opt$exp_file
SAMPLES_FILE <- opt$samples_file
TUMOR_DIR_MAIN <- opt$tumor_dir
OUTPUT_FILE_EXPRESSION <- opt$output_file_expression
OUTPUT_FILE_INDEGREE <- opt$output_file_indegree


source("workflow/bin/tsne_fn.R")

# load expression data
samples <- fread(SAMPLES_FILE)
exp_all <- fread(EXPRESSION_FILE)
colnames(exp_all)[-1] <- samples$sample_id
exp_all <- exp_all[,-1]

# load and combine indegree data
ind_all <- combine_indegree(TUMOR_DIR_MAIN) 
ind_all <- ind_all[,-1]

# run TSNE
tsne_res_ind <- runTSNE_withPCA(ind_all, perplexity = 20, n_pcs = 50)
tsne_res_ind$tsne <- "indegree"
tsne_res_ind$cancer <- samples$cancer[match(tsne_res_ind$id, samples$sample_id)]
write.table(tsne_res_ind,
        file = OUTPUT_FILE_INDEGREE,
        sep = "\t",
        row.names = F,
        col.names = T,
        quote = F)

tsne_res_exp <- runTSNE_withPCA(exp_all, perplexity = 20, n_pcs = 50)
tsne_res_exp$tsne <- "expression"
tsne_res_exp$cancer <- samples$cancer[match(tsne_res_exp$id, samples$sample_id)]
write.table(tsne_res_exp,
        file = OUTPUT_FILE_EXPRESSION,
        sep = "\t",
        row.names = F,
        col.names = T,
        quote = F)
