required_libraries <- c("data.table", "ggplot2")
for (lib in required_libraries) {
  suppressPackageStartupMessages(library(lib, character.only = TRUE,
                                quietly = TRUE))
}

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
