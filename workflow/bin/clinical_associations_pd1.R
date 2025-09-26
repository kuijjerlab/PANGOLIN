#####################
## Load R Packages ##
#####################
required_libraries <- c(
    "data.table", "dplyr", "doParallel", "survminer", "survival",
    "magrittr", "gridExtra", "ggplot2", "ggrepel", "smplot2",
    "effectsize", "rstatix", "optparse"
)

for (lib in required_libraries) {
    suppressPackageStartupMessages(
        library(lib, character.only = TRUE, quietly = TRUE)
    )
}

####################
## Parse Arguments ##
####################
option_list <- list(
    make_option(c("-c", "--cox_results_file"), 
        type = "character",
        default = NULL, help = "Path to the cox results for the PD1 pathway",
        metavar = "character"),
    make_option(c("-g", "--p_threshold"), 
        type = "numeric",
        default = 0.05, help = "P-value threshold for significance",
        metavar = "numeric"),
    make_option(c("-p", "--pd1_scores_file"),
        type = "character",
        default = NULL, help = "Path to the PD1 scores file",
        metavar = "character"),
    make_option(c("-l", "--tumor_clinical_file"), 
        type = "character",
        default = NULL, help = "Path to the tumor clinical file",
        metavar = "character"),
    make_option(c("-n", "--cancer"),
        type = "character",
        default = NULL, help = "Cancer type to process",
        metavar = "character"),
    make_option(c("-r", "--results_pd1_groups"), 
        type = "character",
        default = NULL, help = "Path to the result file 
            comparing clinical groups",
        metavar = "character"),
    make_option(c("-s", "--results_pd1_numeric"),
        type = "character",
        default = NULL, help = "Path to the result file 
            with numerical features",
        metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

COX_RESULTS_FILE <- opt$cox_results_file
P_THRESHOLD <-  opt$p_threshold
PD1_SCORES_FILE <- opt$pd1_scores_file
TUMOR_CLINICAL_FILE <- opt$tumor_clinical_file
CANCER <- opt$cancer
RESULTS_PD1_GROUPS <- opt$results_pd1_groups
RESULTS_PD1_NUMERIC <- opt$results_pd1_numeric

# Validate required arguments
if (is.null(COX_RESULTS_FILE) ||
    is.null(PD1_SCORES_FILE) ||
    is.null(TUMOR_CLINICAL_FILE) ||
    is.null(CANCER) || 
    is.null(RESULTS_PD1_GROUPS) ||
    is.null(RESULTS_PD1_NUMERIC)) {
    stop("Error: All required arguments must be provided.")
}

########################
## Load Helper Scripts ##
########################
source("workflow/bin/clinical_association_pd1_fn.R")
source("workflow/bin/cox_regression_tumor_fn.R")
source("workflow/bin/extract_clinical_data_fn.R")
set.seed(1234)

df <- fread(COX_RESULTS_FILE)

# Hardcode PFI cancers to avoid parsing issues
pfi_cancer <- c("BRCA", "LGG", "PRAD", "READ", "TGCT", "THCA", "THYM")

cat("Using p-value threshold:", P_THRESHOLD, "\n")


df <- df %>%
    filter(
        (type == "PFI" & cancer %in% pfi_cancer) |
        (type == "OS" & !cancer %in% pfi_cancer)
    )

df_clean <- df[df$pvalue < P_THRESHOLD,]

current_cancer_data <- df_clean[df_clean$cancer == toupper(CANCER), ]

if (nrow(current_cancer_data) == 0) {
    cat(
        "Cancer", CANCER,
        "doesn't meet criteria or has no significant results.",
        "Creating empty files.\n"
    )
    writeLines(
        paste(
            "group1", "group2", "p_value", "effect_size", "clin_feature",
            "padjust", "log10padjust", "cancer", "principal_component",
            sep = "\t"
        ),
        RESULTS_PD1_GROUPS
    )
    writeLines(
        paste(
            "clinical_feature", "cor", "pval", "component", 
            "padjust", "log10padjust", "cancer",
            sep = "\t"
        ),
        RESULTS_PD1_NUMERIC
    )
} else {
    cat(
        "Processing", nrow(current_cancer_data),
        "significant components for cancer", CANCER, "\n"
    )
    res_groups_all <- NULL
    res_numeric_all <- NULL

    for (i in seq_len(nrow(current_cancer_data))) {
        component <- current_cancer_data$component[i]

        res_groups <- clin_association_groups_pd1(
            tumor = CANCER,
            clin_cancer_file = TUMOR_CLINICAL_FILE,
            pd1_scores_file = PD1_SCORES_FILE,
            component = component
        )
        res_groups$principal_component <- component
        res_groups$cancer <- CANCER

        res_numeric <- clin_association_numeric_pd1(
            tumor = CANCER,
            clin_cancer_file = TUMOR_CLINICAL_FILE,
            pd1_scores_file = PD1_SCORES_FILE,
            component = component,
            correlation_type = "spearman"
        )
        res_numeric$principal_component <- component
        res_numeric$cancer <- CANCER

        res_groups_all <- rbind(res_groups_all, res_groups)
        res_numeric_all <- rbind(res_numeric_all, res_numeric)
    }

    write.table(
        res_groups_all, RESULTS_PD1_GROUPS, col.names = TRUE,
        row.names = FALSE, sep = "\t", quote = FALSE
    )
    write.table(
        res_numeric_all, RESULTS_PD1_NUMERIC, col.names = TRUE,
        row.names = FALSE, sep = "\t", quote = FALSE
    )
}
