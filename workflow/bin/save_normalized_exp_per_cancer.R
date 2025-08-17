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
        c("-o", "--output_dir"),
        type = "character",
        default = NULL,
        help = "Path to the output directory to save expression files.",
        metavar = "character"
    )
)
opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)


EXPRESSION_FILE <- opt$expression_file
GROUP_FILE <- opt$group_file
OUTPUT_DIR <- opt$output_dir

source("bin/analyze_batch_fn.R")

cat("Loading expression file (this may take a while)...\n")
log2exp <- log2_transform_snail(EXPRESSION_FILE)
cat(sprintf(
    "Loaded expression matrix: %d genes x %d samples\n",
    nrow(log2exp), ncol(log2exp)
))
cat("Loading group file...\n")
groups <- fread(GROUP_FILE, head = FALSE)
cancer_types <- unique(groups$V2)
cat("=== Processing all cancer types ===\n")
cat(sprintf(
    "Cancer types to process: %s\n",
    paste(cancer_types, collapse = ", ")
))

results <- list()
for (cancer in cancer_types) {
    cat(sprintf("\n--- Processing %s ---\n", cancer))
    tryCatch({
        cancer_samples <- groups[groups$V2 == cancer, ]
        if (nrow(cancer_samples) == 0) {
            cat(sprintf(
                "WARNING: No samples found for %s, skipping...\n", cancer
            ))
            results[[cancer]] <- "no_samples"
            next
        }
        cat(sprintf(
            "Found %d samples for %s\n", nrow(cancer_samples), cancer
        ))
        exp_proj <- log2exp[, colnames(log2exp) %in% cancer_samples$V1, drop = FALSE]
        cat(sprintf(
            "Final matrix: %d genes x %d samples\n",
            nrow(exp_proj), ncol(exp_proj)
        ))
        output_file <- file.path(
            OUTPUT_DIR, 
            paste0("normalized_expression_", cancer, ".RData")
        )
        dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
        cat(sprintf("Saving to: %s\n", output_file))
        save(exp_proj, file = output_file)
        results[[cancer]] <- "success"
        cat(sprintf("✓ Successfully processed %s\n", cancer))
    }, error = function(e) {
        cat(sprintf(
            "✗ Error processing %s: %s\n", cancer, e$message
        ))
        results[[cancer]] <- paste("error:", e$message)
    })
}
cat("\n=== Processing Summary ===\n")
for (cancer in cancer_types) {
    cat(sprintf("%s: %s\n", cancer, results[[cancer]]))
}
cat("=== All cancer types processed ===\n")
