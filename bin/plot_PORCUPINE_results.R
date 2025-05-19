library(ComplexHeatmap)
library(data.table)
library(ggplot2)
library(circlize)
library(cowplot)
library(grid)
library(plyr)
library(dplyr)

#####################
## Load R Packages ##
#####################
# List of required libraries
required_libraries <- c("data.table", "optparse", "ComplexHeatmap", "ggplot2",
                        "circlize", "cowplot", "grid", "plyr", "dplyr")
# Load each library and suppress startup messages
for (lib in required_libraries) {
  suppressPackageStartupMessages(library(lib, character.only = TRUE,
                                quietly = TRUE))
}

# Set a random seed for reproducibility
set.seed(1234)


tcga_dir <- "/storage/kuijjerarea/tatiana/DOBERMAN/"
res_dir <- file.path(tcga_dir, "/pan_cancer_results/")

int_res_dir <- file.path(res_dir, "intermediate")
scripts_dir <- file.path(res_dir, "scripts")
fig_dir <- file.path(res_dir, "figs")
pcp_dir <- file.path(res_dir, "PORCUPINE_results_norm")
pcp_dir_mod <- file.path(res_dir, "PORCUPINE_results_modified")
pcp_files <- list.files(pcp_dir, pattern = "pcp_results")
cancers <- gsub("pcp_results_|.txt", "", pcp_files)
msigdb_dir <- file.path(res_dir, "MsigDB")
source(file.path(scripts_dir, "analyze_pcp_results_fn.R"))


# load the filtered set of pathways, and biological functions
res_all <- fread(file.path(int_res_dir, "filtered_PORCUPINE_results.txt"))
ptws_dat <- fread(file.path(int_res_dir,
                "functional_categories_filtered_PORCUPINE_pathways.txt"))
n_pathways <- table(res_all$cancer)
res_list <- split(res_all$pathway, res_all$cancer)
intersection_matrix <- intersect_pathways(res_list, type = c("jaccard"))
colnames(intersection_matrix) <- toupper(colnames(intersection_matrix))
rownames(intersection_matrix) <- toupper(rownames(intersection_matrix))

# plot the results
col_fun = colorRamp2(c(0, 1), c("white", "darkblue"))
pathways_counts <- n_pathways 
column_barplot <- columnAnnotation(
            n_pathways = anno_barplot(as.numeric(pathways_counts),
            gp = gpar(fill = "lightgrey"))
)

pdf(file.path(fig_dir, "heatmap_PORCUPINE_pathways_intersection.pdf"),
                width = 10, height = 10)
ht1 <- Heatmap(intersection_matrix, col = col_fun, 
        width = ncol(intersection_matrix) * unit(5, "mm"),
        height = nrow(intersection_matrix) * unit(5, "mm"),
        heatmap_legend_param = list(title = "jaccard \n similarity"),
        column_names_gp = gpar(fontsize = 10),
        row_names_gp = gpar(fontsize = 10),
        show_row_dend = FALSE, show_column_dend = FALSE,
        top_annotation = column_barplot)
print(ht1)
dev.off()

counts <- table(res_all$cancer)
counts_dat <- data.table(counts)
colnames(counts_dat) <- c("cancer", "n_pathways")
head(counts_dat)

res_all$es <- ifelse(res_all$es > 0, 1, 0)
ptws_summary <- dcast(res_all, pathway ~ cancer, value.var = "es")
ptws_summary[is.na(ptws_summary)] <- 0
sums_counts_ptw <- rowSums(ptws_summary[,-1])
ptws_shared <- ptws_summary[sums_counts_ptw > 5, ]
ids <- ptws_summary$pathway[sums_counts_ptw > 5] # pathways that are shared 

ptws_summary$category <- 
                ptws_dat$category[match(ptws_summary$pathway, ptws_dat$pathway)]

ptws_summary_melt <- melt(ptws_summary[,-1])
stat <- ddply(ptws_summary_melt, .(variable, category), function(x) nrow(x))
stat_category <- 
    ddply(ptws_summary_melt, .(variable, category), function(x) sum(x$value))

stat$n_pathways_from_category <- stat_category$V1
stat <- stat[stat$n_pathways_from_category >0, ]

cols <- c(
        "#db4b83", "#5eb847", "#a95ccf", "#b5b32e", "#5e6ed9", "#d9943c",
        "#6099d5", "#c85529", "#4cc2bc", "#cf444b", "#5bbe7f", "#c74ba8",
        "#557e2e", "#d28dcc", "#a8ad5b", "#7760a5", "#896b2c", "#a34c6c",
        "#39855f", "#d6836d", "#8975ca", "#71a659", "#cb5582", "#c5783e",
        "#b3669e", "#98984d")

length(cols)

categories <- unique(stat$category)
col_categories <- rev(cols)[1:length(categories)]
names(col_categories) <- categories

fixed_order <- ddply(stat, .(category), function(x) nrow(x))
fixed_order <- fixed_order[order(fixed_order$V1, decreasing = F),]
stat$category <- factor(stat$category, levels = fixed_order$category)


ptws_shared$category <- 
                ptws_dat$category[match(ptws_shared$pathway, ptws_dat$pathway)]
categories_shared <- unique(ptws_shared$category)
categories_shared <- 
            rev(fixed_order$category[fixed_order$category %in% categories_shared])
ht_list <- list()


# Constants
ROW_NAME_MAX_LENGTH <- 35
DEFAULT_WIDTH_MM <- 6
HEIGHT_MM_SINGLE <- 0.5
HEIGHT_MM_MULTIPLE <- 2
FONT_SIZE_ROW <- 9
FONT_SIZE_COLUMN <- 10
RECT_COLOR <- "grey"
RECT_LWD <- 0.5


# Initialize heatmap list
ht_list <- list()

# Loop through shared categories
for (i in 1:length(categories_shared)) {
        category_pathway <- categories_shared[i]
        ptws_toplot <- ptws_shared %>%
                    dplyr::filter(category %in% category_pathway)
        pathways <- ptws_toplot$pathway %>%
                    gsub("REACTOME_", "", .)
        pathways <- 
                    substr(pathways, 1, ROW_NAME_MAX_LENGTH)
        ptws_toplot <- ptws_toplot %>%
                    dplyr::select(-pathway, -category) %>%
                    as.matrix()
        colnames(ptws_toplot) <- toupper(colnames(ptws_toplot))
        rownames(ptws_toplot) <- pathways
        # Generate color palette for the heatmap
        col_to_use <- col_categories[category_pathway]
        my_palette <- colorRampPalette(c("white", col_to_use))(n = 2)

        # Determine heatmap dimensions based on the number of rows
        heatmap_height_mm <- ifelse(nrow(ptws_toplot) == 1,
                                    HEIGHT_MM_SINGLE,
                                    HEIGHT_MM_MULTIPLE)

        # Create the heatmap
        ht <- Heatmap(
                ptws_toplot,
                col = my_palette,
                cluster_columns = FALSE,
                show_heatmap_legend = FALSE,
                row_names_gp = gpar(fontsize = FONT_SIZE_ROW),
                column_names_gp = gpar(fontsize = FONT_SIZE_COLUMN),
                rect_gp = gpar(col = RECT_COLOR, lwd = RECT_LWD),
                width = 
                    ncol(intersection_matrix) * unit(DEFAULT_WIDTH_MM, "mm"),
                show_row_dend = FALSE,
                row_names_max_width = unit(15, "cm"),
                row_names_side = "left",
                show_column_dend = TRUE)
        # Add to the heatmap list
        ht_list[[i]] <- ht
}

length(ht_list)
ht_list
ht_list_up <- Reduce(`%v%`, ht_list)
draw(ht_list_up)


p <- ggplot(stat, aes(y = category,  x = toupper(variable), 
            size = n_pathways_from_category, col = category)) +
            geom_point() +
            scale_size(range = c(3, 12)) +
            theme_bw() +
            theme(legend.position =  "none",
            axis.text.x = element_text(angle = 90,
                                        vjust = 0.5, hjust = 1, size = 12),
            axis.text.y = element_text(size = 14)) +
            scale_color_manual(values = col_categories) +
            xlab("") +
            ylab("")

# 1. Draw the heatmap as a grob
heatmap_grob <- grid.grabExpr(draw(ht_list_up,
                padding = unit(c(5, 2, 0, 2), "mm")))

# 2. Extract the ggplot as a grob
gg_grob <- ggplotGrob(p)

# 4. Combine the two grobs using `plot_grid` (from cowplot)
combined_plot <- plot_grid(
    ggdraw() + draw_grob(gg_grob),           # ggplot grob
    ggdraw() + draw_grob(heatmap_grob),     # heatmap grob
    ncol = 1,                               # Arrange vertically
    rel_heights = c(0.5, 0.5)               # Adjust height proportions
)

# 5. Add labels if needed
labeled_plot <- ggdraw(combined_plot) +
                draw_plot_label(label = c("A.", "B."),
                                x = c(0.02, 0.02),
                                y = c(0.98, 0.48),
                                hjust = 0, vjust = 1,
                                size = 15)

# To display the plot
pdf(file.path(fig_dir, "combined_figure_PCP_results_shared_categories.pdf"),
                width = 12, height =12)
print(labeled_plot)
dev.off()

