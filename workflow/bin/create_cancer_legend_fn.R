#' Create a legend for a cancer type plot
#'
#' This function generates a custom legend for a ggplot, mapping `cancer_type`
#' to colors and shapes. The legend is placed at the bottom and displayed in
#' a horizontal direction with formatted text styles.
#'
#' @param data A data frame containing `cancer_type`, `colour`, and `shape`
#'        columns.
#' @return A grob object containing the extracted legend.
#' @import ggplot2
#' @export
create_legend <- function(data) {
    # Define constants
    POINT_SIZE <- 5
    LEGEND_SIZE <- 4
    LEGEND_NROW <- 3
    LEGEND_TITLE_SIZE <- 12
    LEGEND_TEXT_SIZE <- 10

    p <- ggplot(data, aes(
        x = cancer_type, y = 1, color = cancer_type, shape = cancer_type
      )) +
      geom_point(size = POINT_SIZE) +
      scale_color_manual(
        values = setNames(data$colour, data$cancer_type), guide = "legend"
      ) +
      scale_shape_manual(
        values = setNames(data$shape, data$cancer_type), guide = "legend"
      ) +
      guides(
        color = guide_legend(
          title = "Cancer Type", override.aes = list(size = LEGEND_SIZE),
          nrow = LEGEND_NROW
        ),
        shape = guide_legend(
          title = "Cancer Type", override.aes = list(size = LEGEND_SIZE),
          nrow = LEGEND_NROW
        )
      ) +
      theme_minimal() +
      theme(
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.title = 
                element_text(size = LEGEND_TITLE_SIZE, face = "bold"),
        legend.text = element_text(size = LEGEND_TEXT_SIZE),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        plot.margin = margin(10, 50, 10, 50, "pt"), 
      )
    # Extract legend
    legend_grobs <- ggplotGrob(p)$grobs
    legend_index <- which(sapply(legend_grobs, 
                          function(x) x$name) == "guide-box")
    legend <- legend_grobs[[legend_index]]

    return(legend)
}



#' Plot t-SNE for all cancers with custom colors and shapes
#'
#' Generates a t-SNE scatter plot for all cancers, using custom colors and shapes
#' from a color file. Handles missing color/shape assignments and suppresses the
#' legend.
#'
#' @param tsne_res_file Path to t-SNE results file (columns: Dim1, Dim2, cancer).
#' @param cancer_color_file Path to color/shape mapping file (columns: cancer_type,
#'        colour, shape).
#' @param point_size Point size for the plot. Default: 0.7.
#' @param point_alpha Point transparency. Default: 0.7.
#' @param default_color Color for missing cancers. Default: "grey".
#' @param default_shape Shape for missing cancers. Default: 16.
#' @param legend_rows Legend rows (if shown). Default: 3.
#' @param base_font_size Base font size. Default: 12.
#'
#' @return Invisibly returns the ggplot object.
#' @import data.table
#' @import ggplot2
#' @export
plot_TSNE_all_cancers <- function(tsne_res_file, cancer_color_file,
            point_size = 0.7, point_alpha = 0.7, default_color = "grey",
            default_shape = 16, legend_rows = 3, base_font_size = 12) {
          cancer_colors <- fread(cancer_color_file)
          res <- fread(tsne_res_file)
          res$cancer <- toupper(gsub("TCGA-", "", res$cancer))
          res$colour <- cancer_colors$colour[
            match(res$cancer, toupper(cancer_colors$cancer_type))]
          res$shape <- cancer_colors$shape[
            match(res$cancer, toupper(cancer_colors$cancer_type))]
          if (any(is.na(res$colour))) {
            warning("Some cancer types in tsne_res_file do not have colors.")
            res$colour[is.na(res$colour)] <- default_color
          }
          if (any(is.na(res$shape))) {
            warning("Some cancer types in tsne_res_file do not have shapes.")
            res$shape[is.na(res$shape)] <- default_shape
          }
          g1 <- ggplot(res, aes(x = Dim1, 
                      y = Dim2, 
                      color = cancer, 
                      shape = cancer)) +
            geom_point(size = point_size, alpha = point_alpha) +
            scale_colour_manual(values = setNames(cancer_colors$colour,
              toupper(cancer_colors$cancer_type))) +
            scale_shape_manual(values = setNames(cancer_colors$shape,
              toupper(cancer_colors$cancer_type))) +
            labs(x = "t-SNE 1", y = "t-SNE 2") +
            theme_minimal(base_size = base_font_size) +
            theme(
              plot.title = element_text(color = "black", size = base_font_size,
                face = "bold.italic"),
              axis.text = element_text(size = base_font_size - 2),
              axis.title = element_text(size = base_font_size),
              legend.position = "none",
              legend.text = element_text(size = base_font_size - 3),
              legend.title = element_blank(),
              legend.key.size = unit(0.5, "cm")
            ) +
            guides(
              color = guide_legend(nrow = legend_rows, byrow = TRUE),
              shape = guide_legend(nrow = legend_rows, byrow = TRUE)
            )
          print(g1)

}