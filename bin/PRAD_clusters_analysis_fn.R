#' Create Data for Cox Regression in Cancer Analysis
#'
#' Prepares a dataset for Cox  modeling by combining
#' clinical and cluster data, handling covariates, and formatting survival
#' outcome variables.
#'
#' @param tumor_clin_file_path Character. Path to clinical data file.
#' @param tumor_pd1_dir Character. Directory containing tumor PD-1 data.
#' @param covariates Character vector. Covariates to include in Cox model.
#' @param cluster_file Character. Path to clustering result file.
#' @param datatype Character. Type of data to use (default: "clusters").
#' @param type_outcome Character. Survival outcome type (default: "PFI").
#' @param cluster_id Character. Identifier for cluster column (default: "k_4").
#'
#' @return A data frame formatted for Cox regression analysis.

create_data_cox_PRAD <- function(tumor_clin_file_path,
                              tumor_pd1_dir,
                              covariates,
                              cluster_file,
                              datatype = c("clusters"),
                              type_outcome = c("PFI"),
                              cluster_id = "k_4") {
        # Combine clinical and cluster data
        data <- combine_info_for_cancer(tumor_clin_file_path,
                                        tumor_pd1_dir,
                                        cluster_file)

        # Create combined data and Cox regression data
        data_combined <- create_data_combined(data[[datatype]], data$clin)
        data_cox <- create_cox_data(data_combined)

        # Prepare covariates for Cox regression
        processed <- prepare_data_cox_covariates(data_cox, covariates)
        data_cox <- processed$data_cox
        covariates <- processed$covariates
        event_col <- paste0(type_outcome, collapse = "")
        time_col <- paste0(type_outcome, ".time", collapse = "")
        data_cox$event <- as.numeric(data_cox[[event_col]])
        data_cox$ToF_death <- as.numeric(data_cox[[time_col]])
        # Clean the Cox regression data
        data_cox <- clean_data_cox(data_cox)
        return(data_cox = list(data_cox = data_cox,
                                covariates = covariates))
}

#' Plot Cox Proportional Hazards Survival Curves for PRAD Clusters
#' 
#' @param tumor_clin_file_path Character. Path to clinical data file.
#' @param tumor_pd1_dir Character. Directory containing tumor PD-1 data.
#' @param covariates Character vector. Covariates to include in Cox model.
#' @param cluster_file Character. Path to clustering result file.
#' @param datatype Character. Type of data to use (default: "clusters").
#' @param type_outcome Character. Survival outcome type (default: "PFI").
#' @param cluster_id Character. Identifier for cluster column (default: "k_4").
#'
#' @return A Kaplan-Meier survival plot is displayed, 
#' showing the survival distribution by clusters, annotated with group sizes.

plot_PRAD_cox_fit <- function(tumor_clin_file_path,
                              tumor_pd1_dir,
                              covariates,
                              cluster_file,
                              datatype = c("clusters"),
                              type_outcome = c("PFI"),
                              cluster_id = "k_4"){
        data_cox <- create_data_cox_PRAD(tumor_clin_file_path,
                            	tumor_pd1_dir,
                                covariates,
                                cluster_file,
                                datatype = c("clusters"),
                                type_outcome = c("PFI"),
                                cluster_id = "k_4")
        data_cox <- data_cox$data_cox
        feature <- sprintf("`%s`", cluster_id)
        surv_formula <- as.formula(paste("Surv(ToF_death, event) ~", feature))
        group_sizes <- data_cox %>%
                group_by(!!sym(cluster_id)) %>%
                summarise(n = n(), .groups = "drop")

        # 4. Dynamically create label values
        group_vals <- pull(group_sizes, !!sym(cluster_id))
        new_labels <- paste(group_vals, "(n =", group_sizes$n, ")")
        fit <- surv_fit(surv_formula, data = data_cox)
        plot_survival_curve(fit = fit, data_cox, new_labels)

}

#' Plot a Kaplan-Meier Survival Curve
#'
#' Generates a Kaplan-Meier survival plot using `ggsurvplot`, with consistent
#' theming and formatting options for customization.
#'
#' @param fit A `survfit` object from the `survival` package.
#' @param data_cox A data frame containing clinical data used in the survival
#'   fit.
#' @param new_labels Character vector of labels for legend groups.
#' @param pval_coord Numeric vector of length 2 for x and y coordinates of the
#'   p-value annotation. Default is `c(4000, 0.9)`.
#' @param pval_size Numeric value for p-value text size. Default is `7.5`.
#' @param line_size Numeric value for line thickness. Default is `1.2`.
#' @param palette Character string or vector for color palette. Default is
#'   `"Dark2"`.
#' @param font_main_size Numeric value for title font size. Default is `20`.
#' @param font_axis_size Numeric value for axis label font size. Default is
#'   `16`.
#' @param tick_label_size Numeric value for tick and legend text size. Default
#'   is `12`.
#' @param legend_position Position of the legend (e.g., `"bottom"`, `"right"`).
#'   Default is `"bottom"`.
#' @param legend_title Title for the legend. Default is `"Group"`.
#' @param x_label Label for the x-axis. Default is
#'   `"Progression free interval (days)"`.
#' @param y_label Label for the y-axis. Default is
#'   `"Survival probability (%)"`.
#'
#' @return A `ggsurvplot` object representing the survival plot.
#' @export


plot_survival_curve <- function(fit,
                                data_cox,
                                new_labels,
                                pval_coord = c(4000, 0.9),
                                pval_size = 7.5,
                                line_size = 1.2,
                                palette = "Dark2",
                                font_main_size = 20,
                                font_axis_size = 16,
                                tick_label_size = 12,
                                legend_position = "bottom",
                                legend_title = "Group",
                                x_label = "Progression free interval (days)",
                                y_label = "Survival probability (%)") {
    g <- ggsurvplot(
        fit,
        data = data_cox,
        size = line_size,
        palette = palette,
        conf.int = FALSE,
        pval = TRUE,
        pval.size = pval_size,
        pval.coord = pval_coord,
        font.main = c(font_main_size, "bold"),
        font.x = c(font_axis_size, "bold"),
        font.y = c(font_axis_size, "bold"),
        legend = legend_position,
        legend.title = legend_title,
        legend.labs = new_labels,
        xlab = x_label,
        ylab = y_label,
        ggtheme = theme_minimal() +
            theme(
                legend.title = element_blank(),
                axis.line = element_line(),
                legend.text = element_text(size = tick_label_size),
                axis.text.x = element_text(size = tick_label_size),
                axis.text.y = element_text(size = tick_label_size)
            )
    )
    return(g)
}



get_degs_drgs_PRAD <- function(tumor_clin_file_path,
                              tumor_pd1_dir,
                              covariates,
                              cluster_file,
                              datatype = c("clusters"),
                              type_outcome = c("PFI"),
                              cluster_id = "k_4",
                              indegree_file,
                              exp_file,
                              samples_file) {
        data_cox <- create_data_cox_PRAD(tumor_clin_file_path,
                        tumor_pd1_dir,
                        covariates,
                        cluster_file,
                        datatype = datatype,
                        type_outcome = type_outcome,
                        cluster_id = cluster_id)
        res <- perform_limma_PRAD(indegree_file,
                                exp_file,
                                samples_file,
                                data_cox)
        return(res)
}

perform_limma_PRAD <- function(indegree_file, 
                        exp_file,
                        samples_file,
                        data_cox) {
        ind <- load_indegree(ind_file)
        colnames(ind)[-1] <- make_bcr_code(colnames(ind)[-1])
        exp <- load_exp("PRAD", exp_file, samples_file)
        colnames(exp)[-1] <- make_bcr_code(colnames(exp)[-1])
        drgs <- run_limma_PRAD(ind, data_cox)
        degs <- run_limma_PRAD(exp, data_cox)
        all_ranks <- list("indegree" = drgs, "expression" = degs)
        return(all_ranks)
}



run_limma_PRAD <- function(data, 
                        data_cox, 
                        cluster_id = "k_4") {
        # Prepare data
        data <- as_tibble(data)
        genes <- as.character(data[[1]])
        covariates <- data_cox$covariates
        data_cox_df <- data_cox$data
        # Build formula for model matrix
        fml <- as.formula(paste("~ 0 +", cluster_id,
                if (length(covariates) > 0) 
                paste("+", paste(covariates, collapse = " + ")) else ""))
        design_matrix <- model.matrix(fml, data = data_cox_df)

        data_clean <- data[, rownames(design_matrix)]
        data_clean <- as.matrix(data_clean)
        rownames(data_clean) <- genes
        # Extract cluster levels
        cluster_levels <- levels(as.factor(data_cox_df[[cluster_id]]))
        comparisons <- combn(cluster_levels, 2, simplify = FALSE)
        # Build contrast string expressions dynamically
        contrast_strs <- sapply(comparisons, function(cmp) {
                paste0(cluster_id, cmp[1], " - ", cluster_id, cmp[2])
        })
        names(contrast_strs) <- sapply(comparisons, function(cmp) {
                paste0("cl", cmp[1], "_cl", cmp[2])
        })
        # Create contrast matrix
        contrast_matrix <- makeContrasts(contrasts = contrast_strs,
                levels = design_matrix)
        # Run limma
        fit <- lmFit(data_clean, design_matrix)
        fit2 <- contrasts.fit(fit, contrast_matrix)
        fit3 <- eBayes(fit2)
        # Collect results for each contrast
        result_list <- lapply(seq_along(contrast_strs), function(i) {
        res <- topTable(fit3, sort.by = "none", n = Inf, coef = i)
        res$cmp <- names(contrast_strs)[i]
        res$gene <- rownames(res)
        return(res)
        })

        # Combine all results
        res_all <- do.call(rbind, result_list)
        return(res_all)
        }


run_fgsea_PRAD <- function(res_all, 
                        gmt_file, 
                        min_size = 10,
                        max_size = 100,
                        padj_cutoff = 0.01, 
                        es_threshold = 0.5) {

        pt <- fgsea::gmtPathways(gmt_file)
        comparisons <- unique(res_all$cmp)

        fgsea_results_list <- lapply(comparisons, function(comparison) {
        data <- subset(res_all, cmp == comparison)
        ranks <- setNames(data$t, data$gene)

        fgseaRes <- fgsea::fgseaMultilevel(pathways = pt,
                                 stats = ranks,
                                 minSize = min_size,
                                 maxSize = max_size)

        if (nrow(fgseaRes) == 0) return(NULL)

        fgseaRes <- fgseaRes[, 1:7]
        fgseaRes <- subset(fgseaRes,
                        padj < padj_cutoff & abs(ES) >= es_threshold)

        if (nrow(fgseaRes) == 0) return(NULL)

        fgseaRes$cmp <- comparison
        return(fgseaRes)
        })

        fgseaRes_all <- do.call(rbind, Filter(Negate(is.null),
                                                fgsea_results_list))
        return(fgseaRes_all)
        }


plot_fgsea_results <- function(fgseaRes_all,
                        sel_cmps = c("cl_1_cl_2", "cl_1_cl_3", "cl_1_cl_4")) {
        fgseaRes_all$pathway <-
                 gsub("REACTOME_", "", fgseaRes_all$pathway)
        fgseaRes_all$cmp <- gsub("clcluster_", "cl_", fgseaRes_all$cmp)
        fgseaRes_all <- fgseaRes_all[fgseaRes_all$cmp %in% sel_cmps, ]
        fgseaRes_all$pathway <- 
                str_trunc(fgseaRes_all$pathway, width = 70, side = "right")
        fgseaRes_all$log10padj <- -log10(fgseaRes_all$padj)
        g <- fgseaRes_all %>%
            group_by(cmp) %>%
            ungroup() %>%
            mutate(cmp = as.factor(cmp),
                   pathway = reorder_within(pathway, NES, cmp)) %>%
            ggplot(aes(pathway, NES)) +
            geom_col(aes(fill = NES < 0)) +
            coord_flip() +
            scale_x_reordered() +
            labs(x = "", y = "Normalized Enrichment Score") +
            theme_minimal() + theme(legend.position = "none") +
            facet_wrap(~cmp, nrow = 3, scales = "free")
        return(g)
}
