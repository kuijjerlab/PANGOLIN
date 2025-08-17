# List of packages to be loaded
required_libraries <- c("optparse", "ggplot2", "gridExtra",
                        "kableExtra", "data.table")

for (lib in required_libraries) {
  suppressPackageStartupMessages(library(lib, character.only = TRUE,
                                quietly = TRUE))
}

####################
## Read arguments ##
####################
### Options
options(stringsAsFactors = FALSE)

### Command line options
option_list <- list(
    optparse::make_option(
        c("-s", "--summary_table_PD1"),
        type = "character",
        default = NULL,
        help = "Path to the summary PD1 table.",
        metavar = "character"),
    optparse::make_option(
        c("-m", "--output_html_file"),
        type = "character",
        default = NULL,
        help = "Path to the output html file.",
        metavar = "character"))

opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

## Initialize variable ##

SUMMARY_TABLE_PD1 <- opt$summary_table_PD1
HTML_OUTPUT_FILE <- opt$output_html_file

# Load the summary table
data <- fread(SUMMARY_TABLE_PD1)


# Helper function to convert "yes"/"no"/other to symbols
convert_yes_no <- function(x) {
  x <- tolower(trimws(x))
  ifelse(x == "yes", "✓",
         ifelse(x == "no", "-", " "))
}

# Columns to process
cols_to_convert <- c("Heterogeneity", 
                     "Clinical outcome", 
                     "Clinical features", 
                     "Omics clusters", 
                     "Immune Infiltration (CIBERSORTx)")

# Apply the conversion function to each column
data[, (cols_to_convert) := lapply(.SD, convert_yes_no), .SDcols = cols_to_convert]

setnames(data, c("Cancer", "Heterogeneity", 
                 "Clinical outcome", "PDL1 expression", 
                 "Clinical features", "Omics clusters", 
                 "Immune Infiltration (CIBERSORTx)"),
               c("Cancer", "Heterogeneity", 
                 "Clinical<br>outcome", "PDL1<br>expression", 
                 "Clinical<br>features", "Omics<br>clusters", 
                 "Immune Infiltration<br>(CIBERSORTx)"))

# Create the HTML table
table_html <- 
    kable(data, format = "html", escape = FALSE, align = 'lccc') %>%
    kable_styling(bootstrap_options = "striped", 
        full_width = FALSE, position = "left")

# Add meta tag to force UTF-8
meta_tag <- '<meta charset="UTF-8">'

custom_css <- "
<style>
  body {
    margin: 0; /* Remove all body margins */
    padding: 0; /* Remove all padding */
  }

  table {
    table-layout: fixed;
    width: 80%;
    font-family: Arial, 'Segoe UI Symbol', sans-serif;
    text-align: left;
    border-collapse: collapse;
    margin-left: 0; /* Align to the left */
  }

  th, td {
    word-wrap: break-word;
    overflow-wrap: break-word;
    text-align: left;
    padding: 4px 6px; /* Compact cell padding */
    font-size: 13px;  /* Optional: slightly smaller text */
  }

  tr {
    line-height: 1.2; /* Reduce row height */
  }
</style>
"

# Combine meta tag, CSS, and table
table_html_with_css <- paste(meta_tag, custom_css, table_html, sep = "\n")

# Save the HTML file with correct encoding
writeLines(enc2utf8(table_html_with_css),
            HTML_OUTPUT_FILE, useBytes = TRUE)

