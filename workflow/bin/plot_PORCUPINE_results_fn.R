#' Plot a Jaccard Similarity Heatmap with Pathway Counts Annotation
#'
#' @param intersection_matrix Numeric matrix. Jaccard similarity values.
#' @param pathways_counts Numeric vector. Pathway counts for barplot annotation.
#' @param col_fun Function. Color mapping function (default: white to darkblue).
#' @param width_per_col_mm Numeric. Width per column in mm (default: 5).
#' @param height_per_row_mm Numeric. Height per row in mm (default: 5).
#' @param legend_fontsize Integer. Font size for legend (default: 10).
#' @param names_fontsize Integer. Font size for row/column names (default: 10).
#' @param barplot_fill Character. Fill color for barplot (default: "lightgrey").
#'
#' @return A ComplexHeatmap heatmap object.
#' @export
#' 
plot_intersection_heatmap <- function(
                        intersection_matrix,
                        pathways_counts,
                        col_fun = colorRamp2(c(0, 1), c("white", "darkblue")),
                        width_per_col_mm = 5,
                        height_per_row_mm = 5,
                        legend_fontsize = 10,
                        names_fontsize = 10,
                        barplot_fill = "lightgrey"
    ) {
        stopifnot(is.matrix(intersection_matrix))
        if (!is.null(pathways_counts)) {
            stopifnot(length(pathways_counts) == ncol(intersection_matrix))
            column_barplot <- columnAnnotation(
                n_pathways = anno_barplot(
                    as.numeric(pathways_counts),
                    gp = gpar(fill = barplot_fill)
                )
            )
        } else {
            column_barplot <- NULL
        }
        ht1 <- Heatmap(
            intersection_matrix,
            col = col_fun,
            width = ncol(intersection_matrix) * unit(width_per_col_mm, "mm"),
            height = nrow(intersection_matrix) * unit(height_per_row_mm, "mm"),
            heatmap_legend_param = list(
                title = "jaccard \n similarity",
                title_gp = gpar(fontsize = legend_fontsize)
            ),
            column_names_gp = gpar(fontsize = names_fontsize),
            row_names_gp = gpar(fontsize = names_fontsize),
            show_row_dend = FALSE,
            show_column_dend = FALSE,
            top_annotation = column_barplot
        )
        return(ht1)
}



#' Preprocess PORCUPINE Results for Shared Pathway Analysis
#'
#' Summarizes pathway sharing across cancers, annotates with biological
#' categories, and prepares data for downstream plotting.
#'
#' @param res_all Data frame. PORCUPINE results with columns 'pathway',
#'   'cancer', and 'es'.
#' @param ptws_dat Data frame. Pathway annotation table with columns 'pathway'
#'   and 'category'.
#' @param min_shared Integer. Minimum number of cancers a pathway must be
#'   present in to be considered shared (default: 5).
#' @param cols colors for the caterogies
#'
#' @return A list with:
#'   \item{categories_shared}{Character vector of shared categories.}
#'   \item{ptws_shared}{Data frame of shared pathways with category annotation.}
#' @export
preprocess_pathway_results <- function(res_all, 
                                    ptws_dat,
                                    min_shared = 5, 
                                    cols = c(
                                        "#db4b83", "#5eb847", "#a95ccf", 
                                        "#b5b32e", "#5e6ed9", "#d9943c",
                                        "#6099d5", "#c85529", "#4cc2bc", 
                                        "#cf444b", "#5bbe7f", "#c74ba8",
                                        "#557e2e", "#d28dcc", "#a8ad5b", 
                                        "#7760a5", "#896b2c", "#a34c6c",
                                        "#39855f", "#d6836d", "#8975ca", 
                                        "#71a659", "#cb5582", "#c5783e",
                                        "#b3669e", "#98984d")) {
        # Count pathways per cancer
        counts <- table(res_all$cancer)
        counts_dat <- data.table(counts)
        colnames(counts_dat) <- c("cancer", "n_pathways")
        # Binarize enrichment scores
        res_all$es <- ifelse(res_all$es > 0, 1, 0)
        # Summarize pathways across cancers
        ptws_summary <- dcast(res_all, pathway ~ cancer, value.var = "es")
        ptws_summary[is.na(ptws_summary)] <- 0
        sums_counts_ptw <- rowSums(ptws_summary[,-1])
        # Filter for shared pathways
        ptws_shared <- ptws_summary[sums_counts_ptw > min_shared, ]
        shared_ids <- ptws_shared$pathway
        # Annotate with biological category
        ptws_summary$category <- ptws_dat$category[
            match(ptws_summary$pathway, ptws_dat$pathway)
        ]
        # Melt for stats
        ptws_summary_melt <- melt(ptws_summary[,-1])
        stat <- ddply(ptws_summary_melt, .(variable, category), nrow)
        stat_category <- 
            ddply(ptws_summary_melt, .(variable, category), function(x) sum(x$value))
        stat$n_pathways_from_category <- stat_category$V1
        stat <- stat[stat$n_pathways_from_category>0,]
        # Category color and order
        categories <- unique(stat$category)
        col_categories <- rev(cols)[1:length(categories)]
        names(col_categories) <- categories
        fixed_order <- ddply(stat, .(category), nrow)
        fixed_order <- fixed_order[order(fixed_order$V1, decreasing = FALSE), ]
        stat$category <- factor(stat$category, levels = fixed_order$category)
        # Annotate shared pathways
        ptws_shared$category <- ptws_dat$category[
            match(ptws_shared$pathway, ptws_dat$pathway)
        ]
        categories_shared <- unique(ptws_shared$category)
        categories_shared <- rev(
            fixed_order$category[fixed_order$category %in% categories_shared]
        )
        res_list <- list(
            categories_shared = categories_shared,
            ptws_shared = ptws_shared,
            stat = stat,
            col_categories = col_categories
        )
        return(res_list)
}

#' Plot Pathway Category Dot Plot
#'
#' Creates a dot plot showing the number of shared pathways per category and
#' cancer type, with all plot parameters configurable (no magic numbers).
#'
#' @param stat Data frame with columns: category, variable (cancer), and
#'   n_pathways_from_category.
#' @param col_categories Named vector of colors for each category.
#' @param size_range Numeric vector of length 2 for point size range.
#'   (default: c(3, 11))
#' @param axis_text_x_size Integer. Font size for x-axis text (default: 12).
#' @param axis_text_y_size Integer. Font size for y-axis text (default: 14).
#' @param axis_text_x_angle Numeric. Angle for x-axis text (default: 90).
#' @param axis_text_x_vjust Numeric. vjust for x-axis text (default: 0.5).
#' @param axis_text_x_hjust Numeric. hjust for x-axis text (default: 1).
#'
#' @return A ggplot2 object.
#' @export
plot_pathway_category_dotplot <- function(
    stat,
    col_categories,
    size_range = c(3, 10),
    axis_text_x_size = 12,
    axis_text_y_size = 14,
    axis_text_x_angle = 90,
    axis_text_x_vjust = 0.5,
    axis_text_x_hjust = 1
) {
    ggplot(stat, aes(
        y = category,
        x = toupper(variable),
        size = n_pathways_from_category,
        col = category
    )) +
        geom_point() +
        scale_size(range = size_range) +
        theme_bw() +
        theme(
            legend.position = "none",
            axis.text.x = element_text(
                angle = axis_text_x_angle,
                vjust = axis_text_x_vjust,
                hjust = axis_text_x_hjust,
                size = axis_text_x_size
            ),
            axis.text.y = element_text(size = axis_text_y_size)
        ) +
        scale_color_manual(values = col_categories) +
        xlab("") +
        ylab("")
}

#' Make a List of Category Heatmaps
#'
#' @param categories_shared Character vector of categories.
#' @param ptws_shared Data frame of pathways to plot.
#' @param col_categories Named vector of colors for each category.
#' @param intersection_matrix Matrix for width scaling.
#' @param row_name_max_length Integer. Max row name length (default: 35).
#' @param default_width_mm Numeric. Default heatmap width (default: 6).
#' @param height_mm_single Numeric. Height for single-row heatmap (default: 0.5).
#' @param height_mm_multiple Numeric. Height for multi-row heatmap (default: 2).
#' @param font_size_row Integer. Row font size (default: 9).
#' @param font_size_column Integer. Column font size (default: 10).
#' @param rect_color Character. Rectangle color (default: "grey").
#' @param rect_lwd Numeric. Rectangle line width (default: 0.5).
#'
#' @return List of ComplexHeatmap objects.
#' @export
make_category_heatmap_list <- function(
    categories_shared,
    ptws_shared,
    col_categories,
    intersection_matrix,
    row_name_max_length = 35,
    default_width_mm = 6,
    height_mm_single = 0.5,
    height_mm_multiple = 2,
    font_size_row = 9,
    font_size_column = 10,
    rect_color = "grey",
    rect_lwd = 0.5
) {
    ht_list <- list()
    for (i in seq_along(categories_shared)) {
        category_pathway <- categories_shared[i]
        ptws_toplot <- ptws_shared %>%
            dplyr::filter(category %in% category_pathway)
        pathways <- ptws_toplot$pathway %>%
            gsub("REACTOME_", "", .)
        pathways <- substr(pathways, 1, row_name_max_length)
        ptws_toplot <- ptws_toplot %>%
            dplyr::select(-pathway, -category) %>%
            as.matrix()
        colnames(ptws_toplot) <- toupper(colnames(ptws_toplot))
        rownames(ptws_toplot) <- pathways
        col_to_use <- col_categories[category_pathway]
        my_palette <- colorRampPalette(c("white", col_to_use))(n = 2)

        # Check for single-value matrices
        unique_values <- unique(as.vector(ptws_toplot))
        if (length(unique_values) == 1) {
            warning(paste("Skipping heatmap for category:", category_pathway,
                          "because it has only one unique value:", unique_values))
            next
        }

        heatmap_height_mm <- ifelse(
            nrow(ptws_toplot) == 1, height_mm_single, height_mm_multiple
        )
        ht <- Heatmap(
            ptws_toplot,
            col = my_palette,
            cluster_columns = FALSE,
            show_heatmap_legend = FALSE,
            row_names_gp = gpar(fontsize = font_size_row),
            column_names_gp = gpar(fontsize = font_size_column),
            rect_gp = gpar(col = rect_color, lwd = rect_lwd),
            width = ncol(intersection_matrix) * unit(default_width_mm, "mm"),
            show_row_dend = FALSE,
            row_names_max_width = unit(15, "cm"),
            row_names_side = "left",
            show_column_dend = TRUE
        )
        ht_list[[i]] <- ht
    }
    ht_list
}

#' Combine Dotplot and Heatmap Grobs
#'
#' @param dotplot_gg ggplot2 object.
#' @param heatmap_ht_list List of ComplexHeatmap objects.
#' @param rel_heights Numeric vector for relative heights (default: c(0.5, 0.5)).
#' @param labels Character vector for plot labels (default: c("A.", "B.")).
#' @param label_x Numeric vector for label x positions (default: c(0.02, 0.02)).
#' @param label_y Numeric vector for label y positions (default: c(0.98, 0.48)).
#' @param label_size Integer. Label font size (default: 15).
#'
#' @return A cowplot object.
#' @export
combine_dotplot_and_heatmap <- function(
    dotplot_gg,
    heatmap_ht_list,
    rel_heights = c(0.5, 0.5),
    labels = c("A.", "B."),
    label_x = c(0.02, 0.02),
    label_y = c(0.98, 0.48),
    label_size = 15
) {
    heatmap_grob <- grid.grabExpr(draw(heatmap_ht_list,
        padding = unit(c(5, 2, 0, 2), "mm")))
    gg_grob <- ggplotGrob(dotplot_gg)
    combined_plot <- plot_grid(
        ggdraw() + draw_grob(gg_grob),
        ggdraw() + draw_grob(heatmap_grob),
        ncol = 1,
        rel_heights = rel_heights
    )
    labeled_plot <- ggdraw(combined_plot) +
        draw_plot_label(
            label = labels,
            x = label_x,
            y = label_y,
            hjust = 0, vjust = 1,
            size = label_size
        )
    labeled_plot
}