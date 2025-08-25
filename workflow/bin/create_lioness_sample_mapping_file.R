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
        c("-n", "--network_dir"),
        type = "character",
        default = NULL,
        help = "Path to the directory containing LIONESS network TXT files.",
        metavar = "character"
    ),
    optparse::make_option(
        c("-s", "--samples_panda_file"),
        type = "character",
        default = NULL,
        help = "Path to the PANDA samples file 
            (TSV with sample_id and cancer columns).",
        metavar = "character"
    ),
    optparse::make_option(
        c("-o", "--output_file"),
        type = "character",
        default = NULL,
        help = "Path to the output LIONESS sample mapping file.",
        metavar = "character"
    )
)

# Parse command line arguments
opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)


NETWORK_DIR <- opt$network_dir
SAMPLES_PANDA_FILE <- opt$samples_panda_file
OUTPUT_FILE <- opt$output_file


cat("Creating LIONESS sample mapping file...\n")
cat(sprintf("Network directory: %s\n", NETWORK_DIR))
cat(sprintf("PANDA samples file: %s\n", SAMPLES_PANDA_FILE))
cat(sprintf("Output file: %s\n", OUTPUT_FILE))


####################
## Process data   ##
####################

# Find all LIONESS network files
cat("Scanning for LIONESS network files...\n")
filelist <- list.files(NETWORK_DIR, pattern = "lioness.*\\.txt$", recursive = TRUE)

if (length(filelist) == 0) {
    cat("Error: No LIONESS network files found in directory:", NETWORK_DIR, "\n")
    cat("Expected pattern: lioness.*.txt\n")
    quit(status = 1)
}

cat(sprintf("Found %d LIONESS network files\n", length(filelist)))

# Create file information table with proper path extraction
cat("Processing file paths and names...\n")
filelist_dat <- data.table(
    "file" = basename(filelist),
    "file_path" = filelist
)

# Extract network numbers and sort by them
numbers <- as.numeric(regmatches(filelist_dat$file, regexpr("[0-9]+", filelist_dat$file)))

if (any(is.na(numbers))) {
    cat("Warning: Some files don't contain numeric identifiers\n")
    cat("Files without numbers will be placed at the end\n")
    numbers[is.na(numbers)] <- max(numbers, na.rm = TRUE) + seq_along(which(is.na(numbers)))
}

filelist_dat <- filelist_dat[order(numbers), ]
cat(sprintf("Ordered %d files by network number\n", nrow(filelist_dat)))

# Load PANDA samples file
cat("Loading PANDA samples information...\n")
samples <- fread(SAMPLES_PANDA_FILE)

# Validate samples file structure
required_columns <- c("sample_id", "cancer")
missing_columns <- setdiff(required_columns, colnames(samples))
if (length(missing_columns) > 0) {
    cat("Error: Missing required columns in samples file:", paste(missing_columns, collapse = ", "), "\n")
    cat("Required columns: sample_id, cancer\n")
    quit(status = 1)
}

cat(sprintf("Loaded %d samples from PANDA file\n", nrow(samples)))

# Generate expected LIONESS network filenames
samples$lioness_net <- paste0("lioness.", seq_len(nrow(samples)), ".txt")

# Match files to samples
cat("Mapping LIONESS files to TCGA samples...\n")
filelist_dat$tcga_id <-
     samples$sample_id[match(filelist_dat$file, samples$lioness_net)]
filelist_dat$cancer <- 
    samples$cancer[match(filelist_dat$file, samples$lioness_net)]

# Check for unmatched files
unmatched_files <- sum(is.na(filelist_dat$tcga_id))
if (unmatched_files > 0) {
    cat(sprintf("Warning: %d LIONESS files could not be matched to samples\n", unmatched_files))
    cat("These files will have NA values for tcga_id and cancer\n")
}

matched_files <- sum(!is.na(filelist_dat$tcga_id))
cat(sprintf("Successfully matched %d files to TCGA samples\n", matched_files))


####################
## Save output    ##
####################

cat("Writing LIONESS sample mapping file...\n")
write.table(filelist_dat, OUTPUT_FILE,
            col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
