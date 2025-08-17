#########################################################################
## Functions performing clinical associations of PD1 PC componnents
##########################################################################

#' Prepare Data for Clinical Association Analysis
#'
#' Combines PD-1 scores and clinical data, applies tumor-specific cleaning,
#' removes columns with too many missing values, and ensures all columns are
#' suitable for downstream analysis.
#'
#' @param tumor Character. The tumor type.
#' @param clin Data frame. The clinical data.
#' @param pd1_dir Character. Directory path where PD-1 score data is stored.
#' @param na_threshold Numeric. Max allowed fraction of NAs in a column
#'   (default: 0.4, i.e., at least 60% non-NA required).
#'
#' @return Data frame. Cleaned and filtered for clinical association analysis.
#' @export
#' 
prepare_data_for_clin_association <- function(tumor,
                                clin_cancer_file, 
                                tumor_pd1_dir, 
                                na_threshold = 0.4) {
        clin_cancer <- combine_clins(clin_cancer_file)
        pd1_scores <- 
            t(load_pd1_generic(tumor_pd1_dir, type = "pd1_scores"))
        pd1_scores <- pd1_scores[, 1:2, drop = FALSE]
        pd1_scores <- as.data.frame(pd1_scores)
        pd1_scores$bcr_patient_barcode <- rownames(pd1_scores)
        data <- merge(pd1_scores, 
                clin_cancer,
                by = "bcr_patient_barcode", all = TRUE)
        data_filt <- data[!is.na(data$PC1), ]
        if (tumor == "SKCM") data_filt <- clean_skcm_clinical(data_filt)
        if (tumor == "SARC") data_filt <- clean_sarc_clinical(data_filt)
        if (tumor == "CESC") data_filt <- clean_cesc_clinical(data_filt)
        data_filt <- remove_irrelevant_values_columns(data_filt)
        data_filt <- 
            data_filt[, colMeans(is.na(data_filt)) <= na_threshold, drop = FALSE]
        data_filt <- data_filt[
            , vapply(data_filt, function(x) 
                length(unique(na.omit(x))) > 1, logical(1L)),
            drop = FALSE
        ]
        data_filt <- add_cluster_to_name(data_filt)
        data_filt <- replace_commas(data_filt)
        data_filt <- convert_factors_in_dataframe(data_filt)
    return(data_filt)
}

#' Clinical Association with PD-1 PC1-PC2 Scores (Categorical Features)
#'
#' Prepares the data and calculates group-wise associations between clinical
#' features and PD-1 principal components for a given tumor type.
#'
#' @param tumor Character. The tumor type.
#' @param clin Data frame. The clinical data.
#' @param pd1_dir Character. Directory path where PD-1 score data is stored.
#' @param component Character vector. Principal components to analyze
#'   (default: c("PC1", "PC2")).
#' @param na_threshold Numeric. Max allowed fraction of NAs in a column
#'   (default: 0.4).
#'
#' @return Data frame with calculated group-wise associations between clinical
#'   features and PD-1 PC1/PC2 scores.
#' @export
#' @examples
#' clin_association_groups_pd1("prad", clin, "data/pd1/", c("PC1", "PC2"))
clin_association_groups_pd1 <- function(tumor,
                            clin_cancer_file,
                            pd1_dir,
                            component = c("PC1", "PC2"),
                            na_threshold = 0.4) {
    #stopifnot(is.character(tumor), is.data.frame(clin), is.character(pd1_dir))
    data_filt <- prepare_data_for_clin_association(tumor,
        clin_cancer_file, pd1_dir, na_threshold = na_threshold)
    res_clin <- calculate_correlations_groups(data_filt, component)
    res_clin$cancer <- tumor
    return(res_clin)
}

#' Clinical Association with PD-1 PC1-PC2 Scores (Numeric Features)
#'
#' Prepares the data and calculates correlations between numeric clinical
#' features and PD-1 principal components for a given tumor type.
#'
#' @param tumor Character. The tumor type.
#' @param clin Data frame. The clinical data.
#' @param pd1_dir Character. Directory path where PD-1 score data is stored.
#' @param component Character vector. Principal components to analyze
#'   (default: c("PC1", "PC2")).
#' @param correlation_type Character. Correlation method to use, either
#'   "pearson" or "spearman" (default: c("pearson", "spearman")).
#' @param na_threshold Numeric. Max allowed fraction of NAs in a column
#'   (default: 0.4).
#'
#' @return Data frame with calculated correlations between numeric clinical
#'   features and PD-1 PC1/PC2 scores.
#' @export
#' @examples

clin_association_numeric_pd1 <- function(tumor,
                                clin_cancer_file,
                                pd1_dir,
                                component = c("PC1", "PC2"),
                                correlation_type = c("pearson", "spearman"),
                                na_threshold = 0.4) {
    data_filt <- prepare_data_for_clin_association(
        tumor, clin_cancer_file, pd1_dir, na_threshold = na_threshold
    )
    res_clin <- calculate_correlations_numeric(
        data_filt, component, correlation_type
    )
    res_clin$cancer <- tumor
    return(res_clin)
}

#' Calculate Correlations between Clinical Features and PD-1 PC1/PC2 Scores
#'
#' This function calculates the correlations between clinical features and 
#' specified principal components (PC1, PC2) of PD-1 scores.
#'
#' @param data The prepared data containing clinical features and PD-1 scores.
#' @param component The principal components to be analyzed.
#'                  Default is c("PC1", "PC2").
#' @param correlation_type The types of correlation to be calculated.
#'                         Options are "spearman" or "pearson".
#'                         Default is c("spearman", "pearson").
#' @return A dataframe with calculated correlations, p-values,
#'         adjusted p-values, and log-transformed adjusted p-values.
#' @export
calculate_correlations_numeric <- function(data, 
                        component = c("PC1", "PC2"), 
                        correlation_type = c("spearman", "pearson")) {
        component_value <- as.numeric(data[, component, drop = TRUE])
        data <- data[, -c(1:2)]
        data <- data %>% select(where(is.numeric))
        res_all <- NULL
        for (j in 1:length(colnames(data))){
            clin_feature <- colnames(data)[j]
            print(clin_feature)
            clin_feature_value <- 
                    data[, clin_feature, drop = TRUE]
            data_cor <- data.table("component_value" = component_value,
                        "clin_feature_value" = clin_feature_value)
            data_cor <- data_cor[complete.cases(data_cor),]
            pval <- cor.test(data_cor$clin_feature_value, 
                            data_cor$component_value,
                            method = correlation_type)$p.value
            cor <- cor.test(data_cor$clin_feature_value, 
                            data_cor$component_value,
                            method = correlation_type)$estimate
            res <- data.table("clinical_feature" = clin_feature, 
                                "cor" = cor,
                                "pval" = pval,
                                "component" = component)
            res_all <- rbind(res_all, res)
            }
            res_all$padjust <- p.adjust(res_all$pval, method = "fdr")
            res_all$log10padjust <- -log10(res_all$padjust)
            return(res_all)
}

#' Calculate Group-wise Associations between Clinical Features and PD-1 PCs
#'
#' Calculates group-wise statistical associations (e.g., effect size, p-value)
#' between categorical clinical features and specified PD-1 principal components.
#'
#' @param data Data frame. Prepared data with clinical features and PD-1 PCs.
#' @param component Character. Principal component to analyze (default: "PC1").
#'
#' @return Data frame with group-wise test results, effect sizes, p-values,
#'   FDR-adjusted p-values, and -log10 adjusted p-values for each feature.
#' @export
calculate_correlations_groups <- function(data,
                                component = "PC1") {
    component_value <- as.numeric(data[[component]])
    data <- data[, -c(1:2), drop = FALSE]
    data <- data %>% dplyr::select(where(is.character))
    res_all <- NULL
    for (clin_feature in colnames(data)) {
        clin_feature_value <- data[[clin_feature]]
        data_cor <- data.table(
            PC = component_value,
            group = clin_feature_value
        )
        data_cor <- data_cor[complete.cases(data_cor), ]
        results <- perform_pairwise_tests(data_cor, "group", "PC")
        results$clin_feature <- clin_feature
        res_all <- rbind(res_all, results)
    }
    res_all$padjust <- p.adjust(res_all$p_value, method = "fdr")
    res_all$log10padjust <- -log10(res_all$padjust)
    return(res_all)
}

#' Plot Clinical Feature vs. PD-1 Principal Component
#'
#' Plots the distribution or correlation of a clinical feature with a PD-1
#' principal component, supporting both categorical and numeric features.
#' All plot parameters are configurable via function arguments.
#'
#' @param tumor Tumor type (character).
#' @param results_categorical Data frame with categorical association results.
#' @param feature_to_plot Clinical feature to plot (character).
#' @param component Principal component to plot (character).
#' @param clin Clinical data frame.
#' @param pd1_dir Directory for PD-1 data.
#' @param outlier_size Numeric. Size of outlier points (default: 0.5).
#' @param jitter_width Numeric. Jitter width (default: 0.2).
#' @param text_size Numeric. Annotation text size (default: 2).
#' @param axis_text_size Integer. Axis text size (default: 10).
#' @param element_size Integer. Title/label size (default: 10).
#' @param significance_threshold Numeric. FDR threshold (default: 0.01).
#' @param point_size Numeric. Point size for numeric plot (default: 1.5).
#' @param step_increase Numeric. Step increase for geom_signif (default: 0.1).
#' @param box_fill_color Character. Fill color for numeric points (default: "#0f993d").
#' @param point_color Character. Outline color for numeric points (default: "white").
#'
#' @return A ggplot object.
#' @export
plot_clin_feature <- function(
                        tumor,
                        results_categorical,
                        feature_to_plot,
                        component,
                        clin_cancer_file,
                        pd1_dir,
                        outlier_size = 0.5,
                        jitter_width = 0.2,
                        text_size = 2,
                        axis_text_size = 10,
                        element_size = 10,
                        significance_threshold = 0.01,
                        point_size = 1.5,
                        step_increase = 0.1,
                        box_fill_color = "#0f993d",
                        point_color = "white") {
    # Prepare data
    data_filt <- prepare_data_for_clin_association(tumor, 
                clin_cancer_file, pd1_dir)
    data_filt <- data_filt[!is.na(data_filt[[feature_to_plot]]), ]
    feature_value <- data_filt[[feature_to_plot]]
    component_value <- as.numeric(data_filt[[component]])
    basic_theme <- theme(
        plot.title = element_text(size = element_size, face = "plain"),
        axis.title = element_text(size = element_size),
        axis.text = element_text(size = axis_text_size),
        legend.title = element_blank(),
        legend.position = "none"
    )
    data_to_plot <- data.table(
        component_value = component_value,
        feature_value = feature_value
    )
    tumor_label <- toupper(tumor)

    if (!is.numeric(feature_value)) {
        res <- results_categorical %>%
            filter(cancer == tumor) %>%
            filter(clin_feature == feature_to_plot) %>%
            filter(principal_component == component)
        res$padjust <- p_format(res$padjust, accuracy = significance_threshold)
        res$sign <- 
            ifelse(res$padjust < significance_threshold, res$padjust, "NS")
        data_cmps <- res[, c("group1", "group2")]
        cmps <- 
            lapply(1:nrow(data_cmps), function(i) as.character(data_cmps[i, ]))
        g1 <- ggplot(data_to_plot, aes(
                y = as.character(feature_value),
                x = component_value,
                fill = as.character(feature_value)
            )) +
            geom_boxplot(outlier.size = outlier_size) +
            geom_signif(
                comparisons = cmps,
                annotation = res$sign,
                step_increase = step_increase,
                textsize = text_size
            ) +
            geom_jitter(
                position = position_jitter(jitter_width),
                size = outlier_size
            ) +
            xlab(component) +
            theme_bw() +
            ylab("") +
            basic_theme +
            guides(fill = "none") +
            ggtitle(paste0(tumor_label, "_", feature_to_plot))
    } else {
        g1 <- ggplot(data = data_to_plot, aes(
                x = component_value,
                y = as.numeric(feature_value)
            )) +
            geom_point(
                shape = 21,
                fill = box_fill_color,
                color = point_color,
                size = point_size
            ) +
            sm_statCorr(corr_method = "spearman") +
            xlab(component) +
            ylab(feature_to_plot) +
            guides(fill = "none") +
            ggtitle(paste0(tumor_label, "_", feature_to_plot)) +
            basic_theme
    }
    print(g1)
    invisible(g1)
}



#' Add Prefixes to Columns Based on Name Patterns
#'
#' Adds a prefix ("cluster_", "grade_", or "batch_") to values in columns whose
#' names contain the substrings "cluster", "grade", or "batch" (case-insensitive).
#' Only non-NA values are modified.
#'
#' @param df A data frame.
#' @return The input data frame with modified columns.
#' @export
add_cluster_to_name <- function(df) {
    patterns <- c(cluster = "cluster", grade = "grade", batch = "batch")
    for (prefix in names(patterns)) {
        cols <- grep(patterns[[prefix]], colnames(df), ignore.case = TRUE,
                     value = TRUE)
        for (col in cols) {
            df[[col]] <- ifelse(is.na(df[[col]]), NA,
                                paste0(prefix, "_", df[[col]]))
        }
    }
    return(df)
}

#' Convert factor columns in a data frame to numeric or character
#'
#' This function iterates through the columns of a data frame and converts columns
#' that are factors to numeric if all factor levels can be coerced to integers,
#' otherwise converts them to character. Non-factor columns are left unchanged.
#'
#' @param df A data frame.
#' @return The input data frame with modified factor columns.

#' @export
#' 
convert_factors_in_dataframe <- function(df) {
    for (col in colnames(df)) {
        if (is.factor(df[[col]])) {
        # Check if factor is integer or character
        levels_vec <- levels(df[[col]])
        levels_as_integers <- suppressWarnings(as.integer(levels_vec))
        if (all(is.na(levels_as_integers))) {
            df[[col]] <- as.character(df[[col]])
            df[[col]] <- trimws(df[[col]])
            } else {
                  df[[col]] <- as.numeric(df[[col]])
            }
        }
    }
  return(df)
}


#' Replace Commas with Dots in Factor Levels
#'
#' Replaces all commas with dots in the levels of factor columns in a data frame.
#' Non-factor columns are left unchanged.
#'
#' @param df A data frame.
#' @return The input data frame with commas replaced by dots in factor levels.
#' @export
replace_commas <- function(df) {
    df[] <- lapply(df, function(x) {
        if (is.factor(x)) {
            levels(x) <- gsub(",", ".", levels(x), fixed = TRUE)
            x <- factor(x)
        }
        x
    })
    return(df)
}


#' Perform Pairwise Mann-Whitney Tests and Calculate Effect Sizes
#'
#' Performs pairwise Wilcoxon (Mann-Whitney) tests between groups and computes
#' rank-biserial effect sizes for each pair. Returns a data frame with p-values
#' and effect sizes for each group comparison.
#'
#' @param data Data frame with grouping and value columns.
#' @param group_col Character. Name of the grouping column.
#' @param value_col Character. Name of the value column.
#'
#' @return Data frame with group comparisons, p-values, and effect sizes.
#' @export
perform_pairwise_tests <- function(data,
                        group_col,
                        value_col) {
    test_formula <- as.formula(paste(value_col, "~", group_col))
    # Pairwise Wilcoxon tests (no p-value adjustment)
    pairwise_results <- data %>%
        pairwise_wilcox_test(
            test_formula,
            p.adjust.method = "none"
        ) %>%
        select(group1, group2, p) %>%
        mutate(p_value = p) %>%
        select(-p)
    # Rank-biserial effect sizes
    effect_sizes <- data %>%
        wilcox_effsize(test_formula) %>%
        select(group1, group2, effsize) %>%
        mutate(effect_size = effsize) %>%
        select(-effsize)
    # Merge results
    results <- left_join(pairwise_results, effect_sizes,
                         by = c("group1", "group2"))
    return(results)
}

# set of functions to clean messy clinical data
clean_skcm_clinical <- function(data) {
        percent_columns <- c("CC.TT.nTotal.Mut", 
                            "DIPYRIM.C.T.nTotal.Mut", 
                            "DIPYRIM.C.T.n.C.T..mut")
        for (column in percent_columns) {
            data[[column]] <- gsub("%", "", data[[column]])
    }
     return(data)
}


clean_cesc_clinical <- function(data) {
    data$histological_type <- ifelse(
            data$histological_type == "Cervical Squamous Cell Carcinoma", 
                                            "squamous cell carcinomas",
            ifelse(
            data$histological_type == "Adenosquamous", 
            "adenosquamous",
            "adenocarcinomas"
        )
    )
    return(data)
}

clean_sarc_clinical <- function(data) {
    columns_to_fix <- c("purity", "pathologic_tumor_size", 
                        "Cancer_DNA_fraction", "Genome_doublings")
    for (column in columns_to_fix) {
            data[[column]] <- suppressWarnings(as.numeric(data[[column]]))
        }
    cols_to_remove <- c("histological_type", "short_histo")
    data <- data[, !colnames(data) %in% cols_to_remove]
    data$ploidy <- suppressWarnings(as.numeric(data$ploidy))
    data$Subclona_genome_fraction <- 
        suppressWarnings(as.numeric(data$Subclona_genome_fraction))
    return(data)
}

remove_irrelevant_values_columns <- function(data) {
    invalid_values <- c("[Not Evaluated]", "NA", "[Unknown]",
                            "Invalid Number", "[Not Available]", 
                            "[Discrepancy]", "n/a", "-")
    for (value in invalid_values) {
        data[data == value] <- NA
    }
    cols_to_remove <- c("bcr_patient_barcode", "birth_days_to",
                        "OS", "OS.time", "PFI",
                        "type", "PFI.time", "PFI",
                        "patient", "patient_barcode", 
                        "sample", "days_to_birth",
                        "days_to_last_followup", "Country", "Batch", 
                        "Tissue.source.site") 
    data <- data[, !colnames(data) %in% cols_to_remove]
    return(data)
}




#' Combine multiple plots into a grid
#'
#' This function takes a list of plots and arranges them into a grid layout
#' across multiple pages if necessary.
#'
#' @param plot_list A list of plot objects to be arranged in the grid.
#' @param nrow Number of rows per page in the grid layout.
#' @param ncol Number of columns per page in the grid layout.
#'
#' @return A list of plot pages, each containing a grid arrangement of the plots.
#' @export
#'
#' 
combine_plots <- function(plot_list, nrow, ncol) {
      num_plots <- length(plot_list)
      num_pages <- ceiling(num_plots / (nrow * ncol))
      plot_pages <- lapply(1:num_pages, function(page) {
      start_idx <- (page - 1) * (nrow * ncol) + 1
      end_idx <- min(page * (nrow * ncol), num_plots)
      do.call(grid.arrange,
            c(plot_list[start_idx:end_idx], nrow = nrow, ncol = ncol))
})
  return(plot_pages)
}


#' Split a string into cancer and feature components

#' @param string A character string with components separated by underscores.
#'
#' @return A data.table with two columns: "cancer" and "feature". 
#' @export
#' 
split_string_function <- function(string) {
      split_string <- strsplit(string, "_")[[1]]
      res <- data.table("cancer" = split_string[1], 
            "feature" = paste(split_string[-1], collapse = "_"))
      return(res)
}

