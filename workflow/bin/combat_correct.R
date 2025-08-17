#####################
## Load R packages ##
#####################
required_libraries <- c("sva", "data.table", "optparse")
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
        c("-b", "--batch_file"),
        type = "character",
        default = NULL,
        help = "Path to the batch information file.",
        metavar = "character"
    ),
    optparse::make_option(
        c("-c", "--clin_file"),
        type = "character",
        default = NULL,
        help = "Path to the clinical file.",
        metavar = "character"),
    optparse::make_option(
        c("-o", "--output_file"),
        type = "character",
        default = NULL,
        help = "Path to the output file.",
        metavar = "character"
    )
)
opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)


EXPRESSION_FILE <- opt$expression_file
GROUP_FILE <- opt$group_file
BATCH_FILE <- opt$batch_file
CLIN_FILE <- opt$clin_file
OUTPUT_FILE <- opt$output_file


source("workflow/bin/analyze_batch_fn.R")

cat("Loading expression file (this may take a while)...\n")

log2exp <- log2_transform_snail(EXPRESSION_FILE)
cat(sprintf(
    "Loaded expression matrix: %d genes x %d samples\n",
    nrow(log2exp), ncol(log2exp)
))
cat("Loading group file...\n")
groups <- fread(GROUP_FILE, head = FALSE)

cat("Loading batch file...\n")
batch_info <- fread(BATCH_FILE, head = TRUE)

clin <- load_clin_rdata(CLIN_FILE)

# UCEC 
ucec_samples <- groups$V1[groups$V2 == c("TCGA-UCEC")]
exp_ucec <- log2exp[, colnames(log2exp) %in% ucec_samples]
exp_combat_ucec <- combat_correct(exp_ucec, batch_info, clin, "platform")
# combat is generating negative values which can be replaced with 0s
exp_combat_ucec[exp_combat_ucec < 0] <- 0

# READ
read_samples <- groups$V1[groups$V2 == c("TCGA-READ")]
exp_read <- log2exp[, colnames(log2exp) %in% read_samples]
exp_combat_read <- combat_correct(exp_read, batch_info, clin, "platform")
# combat is generating negative values which can be replaced with 0s
exp_combat_read[exp_combat_read < 0] <- 0


# COAD
coad_samples <- groups$V1[groups$V2 == c("TCGA-COAD")]
exp_coad <- log2exp[, colnames(log2exp) %in% coad_samples]
exp_combat_coad <- combat_correct(exp_coad, batch_info, clin, "platform")
# combat is generating negative values which can be replaced with 0s
exp_combat_coad[exp_combat_coad < 0] <- 0


# STAD 
stad_samples <- groups$V1[groups$V2 == c("TCGA-STAD")]
batch_stad <- batch_info[batch_info$Samples %in% stad_samples,]
batch_stad$Year <- ifelse(batch_stad$Year == 2011, "A_batch", "B_batch")
exp_stad <- log2exp[, colnames(log2exp) %in% stad_samples]
exp_combat_stad <- combat_correct(exp_stad, batch_stad, clin, "Year")
# combat is generating negative values which can be replaced with 0s
exp_combat_stad[exp_combat_stad < 0] <- 0

# DLBC 
dlbc_samples <- groups$V1[groups$V2 == c("TCGA-DLBC")]
batch_dlbc <- batch_info[batch_info$Samples %in% dlbc_samples,]
batch_dlbc$Year <- ifelse(batch_dlbc$Year == 2012, "A_batch", "B_batch")
exp_dlbc <- log2exp[, colnames(log2exp) %in% dlbc_samples]
exp_combat_dlbc <- combat_correct(exp_dlbc, batch_dlbc, clin, "Year")
# combat is generating negative values which can be replaced with 0s
exp_combat_dlbc[exp_combat_dlbc < 0] <- 0

# replace the original data with combat-modified
samples_to_replace <- c(ucec_samples, read_samples, coad_samples,
                        stad_samples, dlbc_samples)


exp_combat_mod <- cbind(exp_combat_ucec, exp_combat_read, exp_combat_coad ,
                    exp_combat_stad, exp_combat_dlbc)
log2exp_mod  <- log2exp[, !colnames(log2exp) %in% samples_to_replace]
log2exp_update <- cbind(log2exp_mod, exp_combat_mod )
save(log2exp_update, file = OUTPUT_FILE)

