#####################
## Load R packages ##
#####################
required_libraries <- c("data.table", "optparse", "plyr")
for (lib in required_libraries) {
    suppressPackageStartupMessages(
        library(lib, character.only = TRUE, quietly = TRUE)
    )
}

####################
## Read arguments ##
####################
option_list <- list(
    # Input files
    optparse::make_option(
        c("-e", "--expression_file"),
        type = "character",
        default = NULL,
        help = "Path to the normalized batch corrected expression file (RData).",
        metavar = "character"
    ),
    optparse::make_option(
        c("-b", "--batch_file"),
        type = "character",
        default = NULL,
        help = "Path to the batch information file.",
        metavar = "character"
    ),
    optparse::make_option(
        c("-m", "--motif_file"),
        type = "character",
        default = NULL,
        help = "Path to the motif prior file.",
        metavar = "character"
    ),
    optparse::make_option(
        c("-p", "--ppi_file"),
        type = "character",
        default = NULL,
        help = "Path to the PPI prior file.",
        metavar = "character"
    ),
    optparse::make_option(
        c("-s", "--samples_file"),
        type = "character",
        default = NULL,
        help = "Path to the clean samples file.",
        metavar = "character"
    ),
    optparse::make_option(
        c("-r", "--feature_file"),
        type = "character",
        default = NULL,
        help = "Path to the features file.",
        metavar = "character"
    ),
    optparse::make_option(
        c("-g", "--groups_file"),
        type = "character",
        default = NULL,
        help = "Path to the groups file (TSV with sample-cancer mapping).",
        metavar = "character"
    ),
    
    # Output files
    optparse::make_option(
        c("-f", "--output_motif_file_filtered"),
        type = "character",
        default = NULL,
        help = "Path to the output motif file (filtered for PANDA).",
        metavar = "character"
    ),
    optparse::make_option(
        c("-k", "--output_ppi_file_filtered"),
        type = "character",
        default = NULL,
        help = "Path to the output PPI file (filtered for PANDA).",
        metavar = "character"
    ),
    optparse::make_option(
        c("-o", "--output_samples_file_filtered"),
        type = "character",
        default = NULL,
        help = "Path to the output samples file (filtered for PANDA).",
        metavar = "character"
    ),
    optparse::make_option(
        c("-z", "--output_samples_file_filtered_with_cancer_type"),
        type = "character",
        default = NULL,
        help = "Path to output samples file with cancer type specified.",
        metavar = "character"
    ),
    optparse::make_option(
        c("-n", "--output_expression_file_filtered"),
        type = "character",
        default = NULL,
        help = "Path to the output expression file (filtered for PANDA).",
        metavar = "character"
    ),
    
    # Optional parameters
    optparse::make_option(
        c("--min_sample_expression"),
        type = "integer",
        default = 20,
        help = "Minimum number of samples a gene must be expressed in.",
        metavar = "integer"
    )
)

# Parse command line arguments
opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

# Source helper functions
source("workflow/bin/analyze_batch_fn.R")


# Assign parsed arguments to variables
EXPRESSION_FILE <- opt$expression_file
BATCH_FILE <- opt$batch_file
MOTIF_FILE <- opt$motif_file
PPI_FILE <- opt$ppi_file
SAMPLES_FILE <- opt$samples_file
FEATURE_FILE <- opt$feature_file
GROUPS_FILE <- opt$groups_file
OUTPUT_MOTIF_FILE <- opt$output_motif_file_filtered
OUTPUT_PPI_FILE <- opt$output_ppi_file_filtered
OUTPUT_EXPRESSION_FILE <- opt$output_expression_file_filtered
OUTPUT_SAMPLES_FILE <- opt$output_samples_file_filtered
OUTPUT_SAMPLES_FILE_CANCER_TYPE <- opt$output_samples_file_filtered_with_cancer_type
MIN_SAMPLE_EXPRESSION <- opt$min_sample_expression

cat("Starting PANDA files data preparation...\n")
cat(sprintf("Minimum sample expression threshold: %d\n", MIN_SAMPLE_EXPRESSION))

##################
## Load data    ##
##################

# Load sample information
cat("Loading sample information...\n")
samples_clean <- fread(SAMPLES_FILE)
cat(sprintf("Loaded %d clean samples\n", nrow(samples_clean)))

# Load groups file
cat("Loading groups information...\n")
groups <- fread(GROUPS_FILE)
colnames(groups) <- c("sample_id", "cancer_type")


# Load features file
cat("Loading feature information...\n")
load(FEATURE_FILE, features_env <- new.env())
features <- features_env[['features']]


# Load batch info
cat("Loading batch information...\n")
batch_info <- load_batch(BATCH_FILE)

batch_info <- batch_info[batch_info$Tissues %in% c("cancer"), ]


# Load expression file
cat("Loading expression data...\n")
load(EXPRESSION_FILE, exp_env <- new.env())
log2exp <- exp_env[['log2exp_new']]
cat(sprintf("Loaded expression matrix: %d genes x %d samples\n", 
            nrow(log2exp), ncol(log2exp)))

########################
## Filter protein     ##
## coding genes       ##
########################

cat("Filtering for protein coding genes...\n")
log2exp_clean <- log2exp[features$gene_type %in% c("protein_coding"), ]
features_clean <- features[features$gene_type %in% c("protein_coding"), ]
head(features_clean)

cat(sprintf("After protein coding filter: %d genes\n", nrow(log2exp_clean)))

# Handle duplicate gene names
cat("Handling duplicate gene names...\n")
duplicates <- features_clean[
    duplicated(features_clean$gene_name) |
    duplicated(features_clean$gene_name, fromLast = TRUE), ]

# Remove PAR_Y genes (pseudoautosomal region Y chromosome duplicates)
par_y <- duplicates[grep("_PAR_Y", duplicates$gene_id), ]


# Find remaining duplicates after PAR_Y removal
remaining_duplicates <- duplicates[-grep("_PAR_Y", duplicates$gene_id), ]
remaining_duplicates <- remaining_duplicates[
    duplicated(remaining_duplicates$gene_name) |
    duplicated(remaining_duplicates$gene_name, fromLast = TRUE), ]

genes_to_exclude <- c(par_y$gene_id, remaining_duplicates$gene_id)


# Apply gene filtering

log2exp_clean <- log2exp_clean[!features_clean$gene_id %in% genes_to_exclude, ]

rownames(log2exp_clean) <-
        features_clean$gene_name[!features_clean$gene_id %in% genes_to_exclude]

dim(log2exp_clean) #19932

cat(sprintf("After duplicate removal: %d genes\n", nrow(log2exp_clean)))

########################
## Filter samples     ##
########################

cat("Filtering samples...\n")
# Keep only samples in the clean samples list
log2exp_clean <- log2exp_clean[, 
    colnames(log2exp_clean) %in% samples_clean$sample_id]


####

log2exp_clean <- log2exp_clean[rowSums(log2exp_clean) > 0, ] #19556 genes, 10283 samples
non_zero <- apply(log2exp_clean, 1, function(i) sum(i > 0))
log2exp_clean <- log2exp_clean[non_zero >= MIN_SAMPLE_EXPRESSION, ]

cat(sprintf("After min expression filter (>=%d samples): %d genes\n", 
            MIN_SAMPLE_EXPRESSION, nrow(log2exp_clean)))

#########################
## Process priors      ##
#########################

# Load and process motif prior
cat("Processing motif prior...\n")
motif <- fread(MOTIF_FILE)
colnames(motif) <- c("tf", "gene", "prior")

motif_genes <- unique(motif$gene)
motif_tfs <- unique(motif$tf)

cat(sprintf("Motif prior: %d TFs, %d genes, %d interactions\n", 
            length(motif_tfs), length(motif_genes), nrow(motif)))

# Filter motif prior to genes/TFs in expression data
motif_clean <- motif[motif$gene %in% rownames(log2exp_clean), ]
motif_clean <- motif_clean[motif_clean$tf %in% rownames(log2exp_clean), ]

cat(sprintf("Filtered motif prior: %d TFs, %d genes, %d interactions\n",
            length(unique(motif_clean$tf)), 
            length(unique(motif_clean$gene)), 
            nrow(motif_clean)))

# Load and process PPI prior
cat("Processing PPI prior...\n")
ppi <- fread(PPI_FILE)
ppis <- unique(c(ppi$V1, ppi$V2))

cat(sprintf("PPI prior: %d unique proteins, %d interactions\n", 
            length(ppis), nrow(ppi)))

# Filter PPI prior to genes in expression data
ppi_clean <- ppi[ppi$V1 %in% rownames(log2exp_clean), ]
ppi_clean <- ppi_clean[ppi_clean$V2 %in% rownames(log2exp_clean), ]

cat(sprintf("Filtered PPI prior: %d unique proteins, %d interactions\n",
            length(unique(c(ppi_clean$V1, ppi_clean$V2))), 
            nrow(ppi_clean)))

#########################
## Final filtering     ##
#########################

# Keep only genes that appear in at least one prior
genes_to_include <- unique(c(motif_clean$tf, motif_clean$gene,
                             ppi_clean$V1, ppi_clean$V2))

log2exp_clean <- log2exp_clean[rownames(log2exp_clean) %in% genes_to_include, ]

cat(sprintf("Final expression matrix: %d genes x %d samples\n", 
            nrow(log2exp_clean), ncol(log2exp_clean)))

#########################
## Write output files ##
#########################

cat("Writing output files...\n")

# Write motif prior
write.table(motif_clean, OUTPUT_MOTIF_FILE,
            col.names = FALSE, sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("Motif prior written to: %s\n", OUTPUT_MOTIF_FILE))

# Write PPI prior
write.table(ppi_clean, OUTPUT_PPI_FILE,
            col.names = FALSE, sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("PPI prior written to: %s\n", OUTPUT_PPI_FILE))

# Write expression data
write.table(log2exp_clean, OUTPUT_EXPRESSION_FILE,
            col.names = FALSE, sep = "\t", row.names = TRUE, quote = FALSE)
cat(sprintf("Expression data written to: %s\n", OUTPUT_EXPRESSION_FILE))

# Write samples list
samples_output <- data.table(sample_id = colnames(log2exp_clean))
write.table(samples_output, OUTPUT_SAMPLES_FILE,
            col.names = FALSE, sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("Samples list written to: %s\n", OUTPUT_SAMPLES_FILE))

# Write samples with cancer types
cancers <- groups[match(colnames(log2exp_clean), groups$sample_id), ]
colnames(cancers) <- c("sample_id", "cancer_type")
write.table(cancers, OUTPUT_SAMPLES_FILE_CANCER_TYPE,
            col.names = TRUE, sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("Samples with cancer types written to: %s\n", 
            OUTPUT_SAMPLES_FILE_CANCER_TYPE))

# Summary statistics
cat("\n=== PANDA DATA PREPARATION SUMMARY ===\n")
cat(sprintf("Final genes: %d\n", nrow(log2exp_clean)))
cat(sprintf("Final samples: %d\n", ncol(log2exp_clean)))
cat(sprintf("Cancer types: %d\n", length(unique(cancers$cancer_type))))
cat(sprintf("Motif interactions: %d\n", nrow(motif_clean)))
cat(sprintf("PPI interactions: %d\n", nrow(ppi_clean)))
cat("Preparation completed successfully!\n")