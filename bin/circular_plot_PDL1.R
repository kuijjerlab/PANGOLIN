# List of packages to be loaded
required_libraries <- c("circlize", "data.table", "optparse",
            "dplyr", "plyr", "ComplexHeatmap")

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
        c("-c", "--cox_univariate_results_pd1_pathway"),
        type = "character",
        default = NULL,
        help = "Path to filtered univariate Cox results on PD1 pathway.",
        metavar = "character"),
    optparse::make_option(
        c("-x", "--cox_results_multivariate"),
        type = "character",
        default = NULL,
        help = "Path to multivariate Cox results performed on edges to PD1.",
        metavar = "character"),
    optparse::make_option(
        c("-p", "--ppi_file"),
        type = "character",
        default = NULL,
        help = "Path to the protein-protein interaction (PPI) network file.",
        metavar = "character"),
    optparse::make_option(
        c("-m", "--motif_file"),
        type = "character",
        default = NULL,
        help = "Path to the transcription factor (TF) motif binding file.",
        metavar = "character"),
    optparse::make_option(
        c("-t", "--threshold"),
        type = "numeric",
        default = NULL,
        help = "Minimum times a transcription factor must be selected.",
        metavar = "numeric"),
    optparse::make_option(
        c("-o", "--output_file"),
        type = "character",
        default = NULL,
        help = "Path to the output file to save results.",
        metavar = "character"))

opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

## Initialize variable ##
COX_RESULTS_UNIVARIATE <- opt$cox_univariate_results_pd1_pathway 
COX_RESULTS_MULTIVARIATE <- opt$cox_results_multivariate
PPI_FILE <- opt$ppi_file
MOTIF_FILE <- opt$motif_file
THRESHOLD <- opt$threshold
OUTPUT_FILE <- opt$output_file

## Filter and process data ##
results_cox_univariate <- fread(COX_RESULTS_UNIVARIATE)
data <- fread(COX_RESULTS_MULTIVARIATE)

tumors <- unique(results_cox_univariate$cancer)
pfi_cancer <-  c("BRCA", "LGG", "PRAD", "READ", "TGCT", "THCA", "THYM")
data <- data[data$cancer %in% tumors, ]
data <- data %>%
            filter(type == "PFI" & cancer %in% pfi_cancer |
                    type == "OS" & !cancer %in% pfi_cancer) 

data$TF <- sapply(strsplit(data$edges, "_"), function(x) x[1])
data <- ddply(data, .(TF, cancer, type), function(x) nrow(x))
colnames(data)[4] <- c("ntimes")

data_clean <- data %>%
              filter(ntimes >= THRESHOLD)

all_tfs <- ddply(data_clean, .(TF), function(x) {
        data.frame(
                cancer = paste(x$cancer, collapse = ","),
                cancer_count = nrow(x),
                ntimes = paste(x$ntimes, collapse = ",")
        )
})

all_tfs <- all_tfs[order(all_tfs$cancer_count, decreasing = T), ]
# Subset rows where cancer_count == 1 and order by 'cancer'
subset_ordered <- all_tfs[all_tfs$cancer_count == 1, ]
subset_ordered <- subset_ordered[order(subset_ordered$cancer), ]

# Subset rows where cancer_count != 1 (remain untouched)
subset_untouched <- all_tfs[all_tfs$cancer_count != 1, ]

# Combine the two subsets
all_tfs <- rbind(subset_untouched, subset_ordered)

# Calculate frequency of occurrences by 'cancer'
freq_cancer <- as.data.frame(table(data_clean$cancer))
names(freq_cancer) <- c("cancer", "count")
freq_cancer <- freq_cancer[order(-freq_cancer$count), ]

# Calculate frequency of occurrences by 'TF' and summarize cancer occurrences
freq_tfs <- as.data.frame(table(data_clean$TF))
names(freq_tfs) <- c("TF", "count")
freq_tfs <- freq_tfs[order(-freq_tfs$count), ]

data_clean <- data_clean[data_clean$TF %in% freq_tfs$TF, ]
data_cast <- dcast(data_clean, TF ~ cancer, value.var = "ntimes")
data_cast_t <- t(data_cast[,-1])
colnames(data_cast_t) <- data_cast$TF
data_cast_t[is.na(data_cast_t)] <- 0

tfs_order <-  freq_tfs$TF
data_cast_t <- data_cast_t[, match(tfs_order, colnames(data_cast_t))]

ri <- rownames(data_cast_t)
labels_for_plotting <- colnames(data_cast_t)

split <- 1:ncol(data_cast_t)

data <- t(data_cast_t)
colnames(data) <- ri


# read in ppi data
ppi <- fread(PPI_FILE)
ppi_clean <- ppi[ppi$V1 %in% rownames(data), ]
ppi_clean <- ppi_clean[ppi_clean$V2 %in% rownames(data), ]
ppi_clean
links <- ppi_clean[,1:2]
colnames(links) <- c("reg1", "reg2")

# read in motif data
motifs <- fread(MOTIF_FILE)
motifs <- motifs[grep("^CD274$", motifs$V2),]

# Define label colors based on rownames and motifs
label_colors <- ifelse(
        rownames(data) %in% motifs$V1, 
        "#638ccc", "black")

warm_grey_red_palette <- colorRampPalette(
        c("#fcbba1", "#fc9272", "#99000d"))(50)


col_fun1 <- colorRamp2(
        seq(0, 100, length.out = 51), 
        c("white", warm_grey_red_palette))


# Define sector labels
sectors <- rownames(data)
dir.create(dirname(OUTPUT_FILE), recursive = TRUE, showWarnings = FALSE)
pdf(OUTPUT_FILE, width = 11, height = 11)
circos.clear()
MARGIN_ADJUSTMENT <- 0.1
DEFAULT_GAP <- 0
FINAL_GAP <- 15
START_DEGREE <- 90
TRACK_HEIGHT <- 0.6
CELL_LINE_WIDTH <- 0.2
DEFAULT_BORDER_COLOR <- "black"
LINK_COLOR <- "#638ccc"
DEFAULT_LINK_COLOR <- "grey"
LABEL_ADJUSTMENT <- 0.3
# LABEL_SIZE <- 0.9
LABEL_SIZE <- 0.5
SAMPLE_LABEL_SIZE <- 0.65
SAMPLE_LABEL_ADJUST <- 0.3
LAST_TRACK_INDEX <- length(sectors)
X_LABEL_OFFSET <- 2

# Set up the margins so labels can fit in ##
circos.par(track.margin = rep(MARGIN_ADJUSTMENT, 2)) 
circos.par(
    gap.after = c(rep(DEFAULT_GAP, LAST_TRACK_INDEX - 1), FINAL_GAP),
    start.degree = START_DEGREE
)

circos.heatmap(
    data,
    split = 1:length(sectors),
    cluster = FALSE,
    col = col_fun1,
    track.height = TRACK_HEIGHT,
    cell.lwd = CELL_LINE_WIDTH,
    cell.border = DEFAULT_BORDER_COLOR
)

# Reset the margins so the links will work out
circos.par(track.margin = rep(-MARGIN_ADJUSTMENT, 2)) 

for (i in 1:nrow(links)) {
    from_index <- which(sectors == links$reg1[i])
    to_index <- which(sectors == links$reg2[i])
    circos.heatmap.link(
        from_index,
        to_index,
        col = ifelse(
            links$reg1[i] %in% motifs$V1 | links$reg2[i] %in% motifs$V1,
            LINK_COLOR,
            DEFAULT_LINK_COLOR
        )
    )
}

circos.track(
    track.index = get.current.track.index(),
    panel.fun = function(x, y) {
        circos.text(
            CELL_META$xcenter,
            CELL_META$cell.ylim[2] + convert_y(1, "mm") * 0.9,
            labels_for_plotting[CELL_META$sector.numeric.index],
            col = label_colors[CELL_META$sector.numeric.index],
            facing = "clockwise",
            cex = LABEL_SIZE,
            adj = c(-LABEL_ADJUSTMENT, 0),
            niceFacing = TRUE
        )
    },
    bg.border = NA
)

# Add sample labels to the last sector
circos.track(
    track.index = get.current.track.index(),
    panel.fun = function(x, y) {
        if (CELL_META$sector.numeric.index == LAST_TRACK_INDEX) {
            cn <- toupper(rev(colnames(data)))
            n <- length(cn)
            circos.text(
                rep(CELL_META$cell.xlim[2], n) + convert_x(X_LABEL_OFFSET, "mm"),
                1:n - 0.5,
                cn,
                cex = SAMPLE_LABEL_SIZE,
                adj = c(0, SAMPLE_LABEL_ADJUST),
                facing = "inside"
            )
        }
    },
    bg.border = NA
)

dev.off()
