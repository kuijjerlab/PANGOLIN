library(data.table)
library(dplyr)
library(stringr)
library(ggplot2)
library(cowplot)
library("smplot2")
setwd("/storage/kuijjera")
setwd("/storage/kuijjerarea/tatiana/PANGOLIN/")
cox_results_file <- "data_all/cox_results_all/PD1_pathway_cox_univariate_model_summary_all.txt"
cancer_dir <-   "data_individual_cancers"
source("bin/merge_patient_data_fn.R")
source("bin/plotting_fn.R")

cox_res <- read_all_coxph_results(cox_results_file,
                                pval_threshold = 0.05)
combined_data <- merge_patient_data_all_cancers(cancer_dir)
combined_data
i=1



ltest <- list.files(dir_ind, pattern = "combined_patient_data", recursive = TRUE, full.names = TRUE) %>%
    lapply(function(file) {
        cancer_type <- str_extract(basename(file), "(?<=combined_patient_data_)[A-Z]+")
        df <- fread(file)
        df$cancer <- cancer_type
        return(df)
    }) %>%
    bind_rows() %>%
    filter(cancer %in% coxph_results$cancer) 

predicted_scores <- ltest
plot_list <- lapply(1:nrow(coxph_results), function(i) {
    tumor <- coxph_results$cancer[i]
    component_type <- coxph_results$component_type[i]
    pc_component <- coxph_results$component[i]
    scores_sel <- predicted_scores %>%
                filter(cancer == tumor)
    selected_columns <- c("CD274_exp", pc_component, component_type, "cancer")
    data_sel <- predicted_scores %>%
                filter(cancer == tumor) %>%
                select(all_of(selected_columns))
    colnames(data_sel) <- c("CD274_exp", pc_component, "risk_score", "cancer")
    data_sel <- as.data.frame(data_sel)

    cor_coef <- cor.test(data_sel$CD274_exp, 
                        data_sel[, pc_component],
                            method = "spearman")$estimate
    if(cor_coef < 0){
            data_sel[, pc_component] <- - data_sel[, pc_component]
            } else {
            data_sel[, pc_component] <- data_sel[, pc_component]
            }
        create_pc_cd274_plot(data_sel, pc_component, tumor)

}
)

# Determine the number of plots per page and calculate number of pages needed
plots_per_page <- 18
num_pages <- ceiling(length(plot_list) / plots_per_page)

# Open PDF device
pdf("figs/pdl1_exp_pc1_pc2_selected_components.pdf", 
                width = 12, height = 12)

# Loop through each page
for (page in 1:num_pages) {
        # Determine the index of the plots for the current page
        plot_indices <- 
                ((page - 1) * plots_per_page + 1):min(page * plots_per_page, length(plot_list))
        # Extract the plots for the current page
        page_plots <- plot_list[plot_indices]
        # Combine plots for this page
        if (length(page_plots) > 0) {
        combined_plot <- plot_grid(plotlist = page_plots, ncol = 5)
        print(combined_plot)
        }
        }
dev.off()

