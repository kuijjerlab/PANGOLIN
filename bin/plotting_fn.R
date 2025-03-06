#' Create a scatter plot of principal component vs. CD274 expression
#'
#' This function generates a scatter plot to visualize the relationship between a 
#' principal component and CD274 (PD-L1) expression, with color indicating risk scores.
#'
#' @param data A data frame containing the necessary columns for plotting.
#' @param pc_component A string specifying the column name of the principal component in `data`.
#' @param cancer A string representing the cancer type, used as the plot title.
#'
#' @return A ggplot2 object representing the scatter plot.
#' @export
create_pc_cd274_plot <- function(data, pc_component, cancer) {
    # Validate required columns
    required_columns <- c(pc_component, "CD274_exp", "risk_score")
    missing_cols <- setdiff(required_columns, names(data))
    if (length(missing_cols) > 0) {
        stop(paste("Error: Missing columns:", paste(missing_cols, collapse = ", ")))
    }

    # Compute median for color scaling
    midpoint_value <- median(na.omit(data$risk_score))

    # Create scatter plot
    ggplot(data = data, aes(x = .data[[pc_component]], y = CD274_exp)) +
        geom_point(aes(color = risk_score), size = 2) +
        scale_color_gradient2(low = "blue", mid = "white", high = "red",
                              midpoint = midpoint_value) +  
        theme_minimal() +
        sm_statCorr(corr_method = "spearman", color = "black") +
        ylab("PDL1 (CD274) expression") +
        ggtitle(cancer)
}

#' Generate scatter plots for principal components vs. CD274 expression
#'
#' This function iterates over the results of a Cox proportional hazards model 
#' and generates scatter plots showing the relationship between a principal 
#' component and CD274 expression for different cancer types.
#'
#' @param coxph_results A data frame containing `cancer`, `component_type`, 
#'   and `component` columns.
#' @param predicted_scores A data frame with `CD274_exp`, `cancer`, and principal 
#'   component scores.
#'
#' @return A list of ggplot2 objects representing the generated scatter plots.
#' @export
generate_pc_cd274_plots <- function(cox_results_file,
                            cancer_dir, 
                            pval_threshold = 0.05) {
    cox_res <- read_all_coxph_results(cox_results_file,
                                pval_threshold = 0.05)
    combined_data <- merge_patient_data_all_cancers(cancer_dir)
    plots <- lapply(1:nrow(cox_res), function(i) {
        tumor <- cox_res$cancer[i]
        component_type <- cox_res$component_type[i]
        pc_component <- cox_res$component[i]
        selected_columns <- 
            c("CD274_exp", pc_component, component_type, "cancer")
        # Filter and rename the required columns
        data_sel <- combined_data %>%
            filter(cancer == tumor) %>%
            select(all_of(selected_columns)) %>%
            rename(risk_score = !!sym(component_type)) %>%
            as.data.frame()
        # Compute correlation coefficient
        cor_result <- cor.test(data_sel$CD274_exp, data_sel[[pc_component]],
                               method = "spearman")
        cor_coef <- cor_result$estimate

        if (cor_coef < 0){
            data_sel[, pc_component] <- - data_sel[, pc_component]
                } else {
            data_sel[, pc_component] <- data_sel[, pc_component]
            }

        # Generate and return the plot
        create_pc_cd274_plot(data_sel, pc_component, tumor)
    })
    # Remove NULL values in case of skipped plots
    plots <- Filter(Negate(is.null), plots)
    return(plots)
}