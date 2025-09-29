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
        c("-d", "--network_dir"),
        type = "character",
        default = NULL,
        help = "Directory containing network .RData files.",
        metavar = "character"
    ),
    optparse::make_option(
        c("-o", "--output_file"),
        type = "character",
        default = NULL,
        help = "Output TSV file for network-cancer mapping.",
        metavar = "character"
    )
)

# Parse command line arguments
opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

NETWORK_DIR <- opt$network_dir
OUTPUT_FILE <- opt$output_file

####################
## Generate mapping ##
####################

cat("Scanning network directory:", NETWORK_DIR, "\n")

# Find all network files
net_files <- list.files(NETWORK_DIR, pattern = "^net_TCGA-.*\\.RData$", 
                        full.names = FALSE)

if (length(net_files) == 0) {
    stop("No network files found matching pattern: net_TCGA-*.RData")
}

cat("Found", length(net_files), "network files\n")

# Extract cancer types from filenames
extract_cancer <- function(filename) {
    # Remove net_TCGA- prefix and .RData suffix
    clean_name <- gsub("^net_TCGA-", "", filename)
    clean_name <- gsub("\\.RData$", "", clean_name)
    # Extract cancer code (remove numeric suffixes)
    cancer <- gsub("[0-9]+$", "", clean_name)
    return(cancer)
}

# Create mapping dataframe
mapping <- data.table(
    net_file = net_files,
    cancer = sapply(net_files, extract_cancer)
)

# Sort by cancer then by filename
mapping <- mapping[order(cancer, net_file)]

cat("Generated mapping for", length(unique(mapping$cancer)), 
    "unique cancer types\n")

# Show summary
cancer_counts <- mapping[, .N, by = cancer][order(cancer)]
cat("Files per cancer type:\n")
print(cancer_counts)

# Write mapping file
cat("Writing mapping to:", OUTPUT_FILE, "\n")
write.table(mapping, OUTPUT_FILE, sep = "\t", row.names = FALSE, 
            quote = FALSE, col.names = TRUE)

cat("Done!\n")