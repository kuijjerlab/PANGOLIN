
library(ggplot2)
library(cowplot)
library(tidyr)
library(data.table)
library(dplyr)

coxph_results <- fread("data_all/cox_results_all/PD1_pathway_cox_univariate_model_summary_all.txt")

# tcga_dir <- "/storage/kuijjerarea/tatiana/DOBERMAN/"
# res_dir <- file.path(tcga_dir, "pan_cancer_results")
# fig_dir <- file.path(res_dir, "figs")
# int_res_dir <- file.path(res_dir, "intermediate")
# panda_primary_dir <- file.path(res_dir, "/panda_input_primary/")
# exp_file <- file.path(panda_primary_dir, "exp_tcga_primary.tsv")
# samples_file <- file.path(panda_primary_dir, "samples_cancers_primary.tsv")
# pd1_dir <- file.path(res_dir, "PD1_analysis_norm")

# source(file.path(res_dir, "/scripts/cox_regression_fn_updated.R"))
source("bin/cox_regression_tumor_fn.R")
exp_file <- file.path("panda_input_primary/exp_tcga_primary.tsv")
samples_file <- file.path("panda_input_primary/samples_cancers_primary.tsv")


pdl1_exp <- extract_gene_expression(exp_file, samples_file, "CD274")
pdl1_exp$cancer <- gsub("TCGA-", "", pdl1_exp$cancer)

#coxph_results <- fread(file.path(int_res_dir, "coxph_model_data_all.txt"))
coxph_results <- coxph_results[coxph_results$pvalue <= 0.05,]

pfi_cancer <-  c("BRCA", "LGG", "PRAD", "READ", "TGCT", "THCA", "THYM")

coxph_results <- coxph_results %>%
            filter(type == "PFI" & cancer %in% pfi_cancer |
                    type == "OS" & !cancer %in% pfi_cancer) 

predicted_scores <- fread("data_all/cox_results_all/PD1_pathway_cox_univariate_predited_risk_scores_all.txt")

plot_list <- lapply(1:nrow(coxph_results), function(i) {
        tumor <- toupper(coxph_results$cancer[i])
        pc_component <- coxph_results$component[i]
        datatype <- coxph_results$type[i]
        pd1_scores <- 
                t(load_pd1_generic(tumor, pd1_dir, type = "pd1_scores")) %>%
                as.data.frame() %>%
                mutate(bcr_patient_barcode = rownames(.),
                CD274_exp = pdl1_exp$CD274[
                        match(bcr_patient_barcode, pdl1_exp$bcr_patient_barcode)
           ])
        pd1_scores$CD274_exp <- as.numeric(pd1_scores$CD274_exp)
        scores_sel <- predicted_scores %>%
                filter(cancer == tolower(tumor),
                        component == pc_component,
                        type == datatype)
        pd1_scores <- pd1_scores %>%
                mutate(risk_score = scores_sel$predicted_risk[
                match(bcr_patient_barcode, scores_sel$bcr_patient_barcode)
                ]) %>%
                filter(!is.na(risk_score))
        pd1_scores$risk_score <- as.numeric(scale(pd1_scores$risk_score))
        cor_coef <- cor.test(pd1_scores$CD274_exp, pd1_scores[, pc_component], 
                              method = "spearman")$estimate
        if(cor_coef < 0){
                pd1_scores[, pc_component] <- - pd1_scores[, pc_component]
                } else {
                pd1_scores[, pc_component] <- pd1_scores[, pc_component]
                }
        create_pc_cd274_plot(pd1_scores, pc_component, tumor)
        }
)



# Determine the number of plots per page and calculate number of pages needed
plots_per_page <- 16  
num_pages <- ceiling(length(plot_list) / plots_per_page)

# Open PDF device
pdf(file.path(fig_dir, "pdl1_exp_pc1_pc2_selected_components.pdf"), 
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
        combined_plot <- plot_grid(plotlist = page_plots, ncol = 4)
        print(combined_plot)
        }
        }
dev.off()


# Create dummy data
data <- data.frame(value = seq(-3, 3, length.out = 100))

# Plot
pdf(file.path(fig_dir, "legend.pdf"), 
    width = 4, height = 4)
p <- ggplot(data, aes(x = value, y = 1, color = value)) +
        geom_point() +  # or geom_tile() depending on your preference
        scale_color_gradient2(
        low = "blue", mid = "white", high = "red", midpoint = 0,
        breaks = c(-3, 0, 3),  # Specify where the labels should appear
        labels = c("Better outcome", "", "Worse outcome")  # Custom labels for each break
        ) + theme(legend.position = "bottom")
p
dev.off()

