libraries <- c("data.table", "dplyr", 
            "ggplot2", "ggrepel", "cowplot")
lapply(libraries, library, character.only = TRUE)
#' Load Cancer Colors
#'
#' Reads a file with cancer types and colors
#'
#' @param cancer_color_file Path to the cancer colors file
#' Must contain `cancer_type` and `colour` columns.  
#' 
#' @return A named color vector where names are uppercase cancer types.  
#' Includes "NS" mapped to "grey" for non-significant cases.  
#'
#' @import data.table
#' @export

load_cancer_colors <- function(cancer_color_file) {
    cancer_colors <- fread(cancer_color_file)
    colors <- setNames(cancer_colors$colour, toupper(cancer_colors$cancer_type))
    return(c(c("NS" = "grey"), colors))
}

#' Load and Filter CoxPH Results
#'
#' Reads a Cox proportional hazards model results file and filters rows
#' based on a specified p-value threshold and cancer type criteria.
#'
#' @param coxtph_results_file Path to the CoxPH results file 
#' Must contain `pvalue`, `type`, `cancer`, and `component` columns. 
#' @param pvalue_th Numeric threshold to filter results (default: 0.05).
#'
#' @return A filtered CoxPH results data table with a new `cancer_component`  
#' column combining  `cancer` and `component` values.  
#'
#' @import data.table
#' @import dplyr
#' @export

load_coxph_results <- function(coxph_results_file, 
                    pvalue_th = 0.05) { 
    if (!file.exists(coxph_results_file)) {
        stop("Error: File does not exist - ", coxph_results_file)
    }
    coxph_results <- tryCatch({
        fread(coxph_results_file)
    }, error = function(e) {
        stop("Error reading file: ", e$message)
    })
    if (!all(c("pvalue", "type", "cancer", "component") 
                        %in% colnames(coxph_results))) {
        stop("Error: Missing required columns in the file.")
    }
    coxph_results <- coxph_results[coxph_results$pvalue < pvalue_th, ]
    pfi_cancer <- c("brca", "lgg", "prad", "read", "tgct", "thca", "thym")
    coxph_results <- coxph_results %>%
        filter((type == "PFI" & cancer %in% pfi_cancer) |
               (type == "OS" & !cancer %in% pfi_cancer))
    coxph_results$cancer_component <- paste(
        toupper(coxph_results$cancer), coxph_results$component, sep = "_"
    )
    return(coxph_results)
}



#' Standardize Column Names for Results Data
#'
#' Renames columns in a results data table to standardized names for  
#' consistency in downstream processing.
#'
#' @param res A data table with six columns to be renamed.
#'
#' @return The input data table with standardized column names.
#'
#' @import data.table
#' @export

setnames_res <- function(res) {
    col_ids <- c("effect_size_cor", "clin_feature", "padjust",
                 "log10padjust", "cancer", "principal_component")
    if (ncol(res) != length(col_ids)) {
        stop("Error: Data table must have exactly ",
             length(col_ids), " columns.")
    }
    setnames(res, colnames(res), col_ids)
    return(res)
}


#' Process Categorical Results
#'
#' Loads and processes categorical results by filtering based on CoxPH  
#' results and selecting the most significant effect size per feature.
#'
#' @param res_categorical_file Path to the categorical results file 
#' @param coxph_results_file Path to the CoxPH results file 
#' @return A processed data table with standardized column names,  
#' filtered by significant effect sizes.
#'
#' @import data.table
#' @import dplyr
#' @export

process_categorical_results <- function(res_categorical_file,
                                        coxph_results_file) {
    if (!file.exists(res_categorical_file)) {
        stop("Error: File does not exist - ", res_categorical_file)
    }
    if (!file.exists(coxph_results_file)) {
        stop("Error: File does not exist - ", coxph_results_file)
    }
    coxph_results <- load_coxph_results(coxph_results_file) 
    res_categorical <- tryCatch({
        fread(res_categorical_file)
    }, error = function(e) {
        stop("Error reading file: ", e$message)
    })

    required_cols <- c("cancer", "principal_component",
                    "clin_feature", "effect_size")
    if (!all(required_cols %in% colnames(res_categorical))) {
        stop("Error: Missing required columns in categorical results file.")
    }
    res_categorical$cancer_component <- paste(
        toupper(res_categorical$cancer), 
        res_categorical$principal_component, sep = "_"
    )
    res_categorical <- res_categorical %>%
        filter(cancer_component %in% coxph_results$cancer_component)
    res_cat <- res_categorical %>%
        mutate(feature_cancer_component = paste(
            cancer, clin_feature, principal_component, sep = "_"
        )) %>%
        group_by(feature_cancer_component) %>%
        filter(effect_size == max(effect_size)) %>%
        slice(1) %>%  # Select first row in case of ties
        ungroup() %>%
        select(effect_size, clin_feature, padjust, 
            log10padjust, cancer, principal_component)
    res_cat <- setnames_res(res_cat)
    return(res_cat)
}

#' Process Numeric Results
#'
#' Loads and processes numeric results by filtering based on CoxPH results  
#' and selecting relevant columns for further analysis.
#'
#' @param res_numeric_file Path to the numeric results file 
#' Must contain `cancer`, `principal_component`, `clinical_feature`, and `cor`.
#' @param coxph_results_file Path to the CoxPH results file
#' Used for filtering numeric results by `cancer_component`.
#'
#' @return A processed data table with standardized column names,  
#' filtered to include only relevant results.
#'
#' @import data.table
#' @import dplyr
#' @export

process_numeric_results <- function(res_numeric_file,
                                    coxph_results_file) {
    if (!file.exists(res_numeric_file)) {
        stop("Error: File does not exist - ", res_numeric_file)
    }
    if (!file.exists(coxph_results_file)) {
        stop("Error: File does not exist - ", coxph_results_file)
    }
    coxph_results <- load_coxph_results(coxph_results_file)
    res_numeric <- tryCatch({
        fread(res_numeric_file)
    }, error = function(e) {
        stop("Error reading file: ", e$message)
    })
    required_cols <- c("cancer", "principal_component", "clinical_feature", "cor")
    if (!all(required_cols %in% colnames(res_numeric))) {
        stop("Error: Missing required columns in numeric results file.")
    }
    res_numeric$cancer_component <- paste(
        toupper(res_numeric$cancer),
        res_numeric$principal_component, sep = "_"
    )
    res_numeric <- res_numeric %>%
        filter(cancer_component %in% coxph_results$cancer_component)
    res_num <- res_numeric %>%
        mutate(feature_cancer_component = paste(
            cancer, clinical_feature, principal_component, sep = "_"
        )) %>%
        select(cor, clinical_feature, padjust, 
        log10padjust, cancer, principal_component)
    res_num <- setnames_res(res_num)
    return(res_num)
}

#' Combine and Process Categorical and Numeric Results
#'
#' Merges categorical and numeric results, standardizes cancer names,  
#' handles infinite values in `log10padjust`, and assigns conditions.
#'
#' @param res_cat A data table of processed categorical results.  
#' @param res_num A data table of processed numeric results.  
#' @param padjust_th Numeric threshold for `padjust` filtering 
#' (default: 0.01).  
#' @param effect_size_th Numeric threshold for `effect_size_cor` filtering  
#' (default: 0.4).
#'
#' @return A merged and processed data table with a new `condition` column  
#' and a `cancer_component` identifier.

combine_data <- function(res_cat, res_num, padjust_th = 0.01,
                         effect_size_th = 0.4) { 
    required_cols <- c("effect_size_cor", "cancer", "principal_component",
                       "log10padjust", "padjust")
    if (!all(required_cols %in% colnames(res_cat))) {
        stop("Error: Missing required columns in categorical results.")
    }
    if (!all(required_cols %in% colnames(res_num))) {
        stop("Error: Missing required columns in numeric results.")
    }
    data_all <- rbind(res_cat, res_num)
    data_all$cancer <- toupper(data_all$cancer)
    max_padjust <- data_all %>%
        filter(!is.infinite(log10padjust)) %>%
        summarise(max_log10padjust = max(log10padjust, na.rm = TRUE)) %>%
        pull(max_log10padjust)
    data_all <- data_all %>%
        mutate(log10padjust = ifelse(
                   is.infinite(log10padjust), max_padjust, log10padjust),
               condition = ifelse(
                   padjust > padjust_th | abs(effect_size_cor) < effect_size_th,
                   "NS", cancer)
        )
    data_all$cancer_component <- paste(
        data_all$cancer, data_all$principal_component, sep = "_"
    )
    return(data_all)
}


#' Create a Volcano Plot for Effect Size and p-Adjusted Values
#'
#' Generates a scatter plot showing the relationship between absolute  
#' effect size and -log10(p-adjusted) values for a given tumor and  
#' principal component, with all plot parameters configurable.
#'
#' @param data A data table containing `effect_size_cor`, `log10padjust`,  
#'   `condition`, and `clin_feature` columns.
#' @param pc_component Principal component identifier for the plot.
#' @param tumor Cancer type used as the plot title.
#' @param colors_to_use A named vector mapping conditions to colors.
#' @param effect_size_threshold Numeric threshold for highlighting features  
#'   (default: 0.4).
#' @param padjust_threshold Numeric threshold for significant p-adjusted  
#'   values (default: 0.01).
#' @param alpha_value Numeric. Point transparency (default: 0.7).
#' @param size_range Numeric vector. Range for point sizes (default: c(0, 4)).
#' @param effect_size_max Numeric. Maximum effect size for scaling (default: 0.85).
#' @param x_limits Numeric vector. X-axis limits (default: c(0, 0.8)).
#' @param y_limits Numeric vector. Y-axis limits (default: c(0, 55)).
#' @param base_font_size Integer. Base font size (default: 10).
#' @param text_size Numeric. Size for text labels (default: 2.5).
#' @param max_overlaps Integer. Max overlaps for text repel (default: 30).
#' @param box_padding Numeric. Padding for text repel box (default: 0.2).
#' @param point_padding Numeric. Padding for text repel point (default: 0.3).
#'
#' @return A ggplot2 object representing the scatter plot.
#' @export
plot_clinical_associations <- function(
                            data,
                            pc_component, 
                            tumor,
                            colors_to_use,
                            effect_size_threshold = 0.4, 
                            padjust_threshold = 0.01,
                            alpha_value = 0.7,
                            size_range = c(0, 4),
                            effect_size_max = 0.85,
                            x_limits = c(0, 0.8),
                            y_limits = c(0, 55),
                            base_font_size = 10,
                            text_size = 2.5,
                            max_overlaps = 30,
                            box_padding = 0.2,
                            point_padding = 0.3
        ) {
            ggplot(data, aes(abs(effect_size_cor),
                            y = log10padjust,
                            color = condition, 
                            fill = condition)) +
                geom_point(aes(size = abs(effect_size_cor)), 
                            alpha = alpha_value) +
                scale_colour_manual(values = colors_to_use) +
                scale_size(range = size_range, limits = c(0, effect_size_max)) +
                scale_x_continuous(limits = x_limits) +
                scale_y_continuous(limits = y_limits) +
                theme_light(base_size = base_font_size) +
                theme(plot.title = element_text(size = base_font_size, 
                                    face = "bold", hjust = 0.5)) +
                geom_text_repel(aes(label = ifelse(
                    abs(effect_size_cor) >= effect_size_threshold & 
                    padjust <= padjust_threshold, clin_feature, "")),
                    size = text_size, max.overlaps = max_overlaps, 
                    box.padding = box_padding, point.padding = point_padding) +
                theme(legend.position = "none",
                    legend.text = element_text(size = base_font_size),
                    axis.text.x = element_text(size = base_font_size),
                    axis.text.y = element_text(size = base_font_size),
                    axis.title.x = element_text(size = base_font_size),
                    axis.title.y = element_text(size = base_font_size),
                    legend.title = element_blank(),
                    panel.border = element_rect(color = "black", 
                                fill = NA),
                    panel.grid.major = element_line(size = 0.5,
                                color = "lightgray"),
                    panel.grid.minor = element_blank(),
                    legend.box.spacing = unit(0.03, "cm"),
                    legend.spacing.x = unit(0.01, "cm")) +
                ylab("-log10(p-adjusted)") + 
                xlab("Effect Size / Correlation") +
                ggtitle(tumor)
}