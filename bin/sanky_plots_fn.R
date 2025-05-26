#' Create Datasets for Sankey Plot and Rand Index Calculation
#'
#' This function processes clustering data for multiple cancer types to create
#' datasets for Sankey plots and calculates the Rand Index (RI) between
#' expression and indegree clusters.
#'
#' @param data A data frame with clustering info. Must include columns:
#'   `cancer`, `bcr_patient_barcode`, `datatype`, and `cluster_id`.
#' @return A list with two elements:
#'   \item{datasets}{List of data frames for Sankey plots, one per cancer type.}
#'   \item{ri_index}{List of data frames with Rand Index (RI) for each cancer.}
#' @export
make_dataset_for_sanky <- function(data) {
    required_cols <- 
        c("cancer", "bcr_patient_barcode", "datatype", "cluster_id")
    if (!all(required_cols %in% colnames(data))) {
        stop("Input must contain: ", paste(required_cols, collapse = ", "))
    }
    cancers <- unique(data$cancer)
    cancers <- cancers[order(cancers)]
    datasets <- list()
    ri_index <- list()
    for (tumor in cancers) {
    data_subset <- data %>%
        filter(cancer == tumor) %>%
        mutate(bcr_code_datatype = 
            paste(bcr_patient_barcode, datatype, sep = "_")) %>%
        distinct(bcr_code_datatype, .keep_all = TRUE) %>%
        mutate(cluster_id = case_when(
            datatype == "indegree" ~ paste0("ind_", cluster_id),
            datatype == "expression" ~ paste0("exp_", cluster_id)
        ))
    data_cast <- data_subset %>%
      reshape2::dcast(bcr_patient_barcode ~ datatype, value.var = "cluster_id")
    data_cast[is.na(data_cast)] <- "not_avail"
    df <- data_cast %>% make_long(expression, indegree)
    RI_exp_ind <- 
        mclust::adjustedRandIndex(data_cast$expression, data_cast$indegree)
    ri_data <- data.table("RI_exp_ind" = RI_exp_ind)
    datasets[[tumor]] <- df
    ri_index[[tumor]] <- ri_data
    }

    clean_datasets <- datasets[sapply(datasets, function(x) !is.null(x))]
    for (name in names(clean_datasets)) clean_datasets[[name]]$cancer <- name

    ri_index <- ri_index[sapply(ri_index, function(x) !is.null(x))]
    for (name in names(ri_index)) ri_index[[name]]$cancer <- name
    return(list(datasets = clean_datasets, ri_index = ri_index))
    }



#' Create a Sankey Diagram for a Given Cancer Dataset
#'
#' Generates a Sankey plot using ggplot2 and ggsankey for a dataset representing
#' node transitions or flows.
#'
#' @param df A data frame containing Sankey plot data. Must include `x`,
#'   `next_x`, `node`, and `next_node` columns.
#' @param cancer Title for the plot (e.g., type of cancer).
#' @param flow_alpha Transparency of flows (0 to 1). Default is 0.5.
#' @param node_color Color or number for node outline. Default is 1 (black).
#' @param label_size Size of node labels. Default is 2.5.
#' @param label_color Text color for labels. Default is 1 (black).
#' @param label_fill Background fill for labels. Default is "white".
#' @param base_size Base font size for theme. Default is 16.
#' @param title_hjust Horizontal justification of title (0 to 1). 
#' Default is 0.5.
#'
#' @return A ggplot object representing the Sankey diagram.
#' @export
sanky_plot_viridis <- function(df,
                               cancer,
                               flow_alpha = 0.5,
                               node_color = 1,
                               label_size = 2.5,
                               label_color = 1,
                               label_fill = "white",
                               base_size = 16,
                               title_hjust = 0.5) {
        g1 <- ggplot(
            df,
            aes(
            x = x,
            next_x = next_x,
            node = node,
            next_node = next_node,
            fill = factor(node),
            label = node
            )
        ) +
            geom_sankey(flow.alpha = flow_alpha, node.color = node_color) +
            geom_sankey_label(
            size = label_size,
            color = label_color,
            fill = label_fill
            ) +
            scale_fill_viridis_d() +
            theme_sankey(base_size = base_size) +
            theme(
            legend.position = "none",
            axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            plot.title = element_text(hjust = title_hjust)
            ) +
            ggtitle(cancer)

        return(g1)
        }
