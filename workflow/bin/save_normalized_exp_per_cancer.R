#####################
## Load R packages ##
#####################
required_libraries <- c("data.table", "optparse")
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
        c("-e", "--expression_file"),
        type = "character",
        default = NULL,
        help = "Path to the normalized expression file.",
        metavar = "character"
    ),
    optparse::make_option(
        c("-g", "--group_file"),
        type = "character",
        default = NULL,
        help = "Path to the sample group file.",
        metavar = "character"
    ),
    optparse::make_option(
        c("-c", "--cancer"),
        type = "character",
        default = NULL,
        help = "Cancer type to process.",
        metavar = "character"
    ),
    optparse::make_option(
        c("-o", "--output_file"),
        type = "character",
        default = NULL,
        help = "Path to the output 
            file for the cancer-specific expression matrix.",
        metavar = "character"
    )
)
opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

EXPRESSION_FILE <- opt$expression_file
GROUP_FILE <- opt$group_file
CANCER <- opt$cancer
OUTPUT_FILE <- opt$output_file

CANCER <- paste0("TCGA-", CANCER)
if (is.null(CANCER) || is.null(OUTPUT_FILE)) {
    stop("Both --cancer and --output_file must be specified.")
}

source("workflow/bin/analyze_batch_fn.R")

cat("Loading expression file (this may take a while)...\n")
log2exp <- log2_transform_snail(EXPRESSION_FILE)
cat(sprintf(
    "Loaded expression matrix: %d genes x %d samples\n",
    nrow(log2exp), ncol(log2exp)
))
cat("Loading group file...\n")
groups <- fread(GROUP_FILE, head = FALSE)

cat(sprintf("\n--- Processing %s ---\n", CANCER))
tryCatch({
    cancer_samples <- groups[groups$V2 == CANCER, ]
    if (nrow(cancer_samples) == 0) {
        stop(sprintf("No samples found for %s, skipping...\n", CANCER))
    }
    cat(sprintf(
        "Found %d samples for %s\n", nrow(cancer_samples), CANCER
    ))
    exp_proj <- log2exp[, 
        colnames(log2exp) %in% cancer_samples$V1, drop = FALSE]
    cat(sprintf(
        "Final matrix: %d genes x %d samples\n",
        nrow(exp_proj), ncol(exp_proj)
    ))
    dir.create(dirname(OUTPUT_FILE), recursive = TRUE, showWarnings = FALSE)
    cat(sprintf("Saving to: %s\n", OUTPUT_FILE))
    save(exp_proj, file = OUTPUT_FILE)
    cat(sprintf("✓ Successfully processed %s\n", CANCER))
}, error = function(e) {
    cat(sprintf("✗ Error processing %s: %s\n", CANCER, e$message))
    stop(e)
})