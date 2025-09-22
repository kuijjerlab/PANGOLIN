## Run univatiate cox model on pd1-pathway based scores for each cancer ##        
rule run_univariate_cox_pd1_pathway:
    input:
        clin_file = TUMOR_CLIN_FILE,
        tumor_pd1_dir = TUMOR_PD1_DIR
    output:
        out_file_summary = OUTPUT_CANCER_UNIVARIATE_COX_SUMMARY,
        out_file_predicted_scores = OUTPUT_CANCER_UNIVARIATE_COX_PREDICTED_SCORES
    log:
        "logs/run_univariate_cox_pd1_pathway_{cancer}.log"
    params:
        bin = config["bin"]
    message:
        "Running univariate Cox model on pd1-pathway based scores for: {wildcards.cancer}"
    shell:
        """
        Rscript {params.bin}/run_univariate_cox_pd1_pathway.R \
            --tumor_clin_file_path {input.clin_file} \
            --tumor_pd1_dir {input.tumor_pd1_dir} \
            --cox_model_summary {output.out_file_summary} \
            --cox_predicted_risk {output.out_file_predicted_scores} \
            > {log} 2>&1
        """
rule combine_pd1_pathway_univariate_summary_results:
    input:
        expand(OUTPUT_CANCER_UNIVARIATE_COX_SUMMARY, cancer=CANCER_TYPES)
    output:
        UNIVARIATE_COX_SUMMARY_ALL
    message:
        "Combining all univariate Cox regression results into one table."
    log:
        "logs/combine_pd1_pathway_univariate_summary_results.log"
    shell:
        """
        echo -e "cancer\\t$(head -n 1 {input[0]})" > {output}
        for file in {input}; do
            cancer=$(basename $(dirname $(dirname $file)))
            tail -n +2 $file | awk -v cancer=$cancer '{{print cancer"\\t"$0}}' >> {output}
        done
        echo "Finished combining results. Output written to: {output}" >> {log}
        """

rule combine_pd1_pathway_univariate_prediction_scores:
    input:
        expand(OUTPUT_CANCER_UNIVARIATE_COX_PREDICTED_SCORES, cancer=CANCER_TYPES)
    output:
        UNIVARIATE_COX_PREDICTED_SCORES_ALL
    log:
        "logs/combine_pd1_pathway_univariate_prediction_scores.log"
    message:
        "Combining all univariate Cox prediction risk results into one table."
    shell:
        """
        echo -e "cancer\\t$(head -n 1 {input[0]})" > {output}
        for file in {input}; do
            cancer=$(basename $(dirname $(dirname $file)))
            tail -n +2 $file | awk -v cancer=$cancer '{{print cancer"\\t"$0}}' >> {output}
        done
        echo "Finished combining results. Output written to: {output}" >> {log}
        """

## Clean the univatiate cox model on pd1-pathway results based of the PORCUPINE results ##        
rule clean_univariate_cox_pd1_pathway:
    input:
        univarite_cox_summary_all = UNIVARIATE_COX_SUMMARY_ALL,
        porcupine_results_all_filtered = PORCUPINE_RESULTS_ALL
    output:
        univarite_cox_summary_all_filtered = UNIVARIATE_COX_SUMMARY_ALL_FILTERED
    log:
        "logs/clean_univariate_cox_pd1_pathway.log"
    params:
        bin = config["bin"]
    message:
        "Cleaning the univariate Cox model on pd1-pathway results based on the PORCUPINE results"
    shell:
        """
        Rscript {params.bin}/clean_pd1_pathway_cox_results_by_porcupine.R \
            --cox_summary_all_cancers {input.univarite_cox_summary_all} \
            --porcupine_filtered_results {input.porcupine_results_all_filtered} \
            --cox_summary_all_cancers_filtered {output.univarite_cox_summary_all_filtered} \
            > {log} 2>&1
        """



### Combine the PD1 pathway heterogeneity scores  with PDL1 expression and PD1 pathway scores ###
### and immune infiltration ###

rule merge_patient_data:
    input:
        tumor_pd1_dir = TUMOR_PD1_DIR,
        risk_score = OUTPUT_CANCER_UNIVARIATE_COX_PREDICTED_SCORES,
        immune_file = IMMUNE_FILE
    output:
        out_file = OUTPUT_COMBINED_PATIENT_DATA_CANCER
    log:
        "logs/merge_patient_data_{cancer}.log"
    message:
        "Merging patient data for: {wildcards.cancer}"
    params:
        bin = config["bin"]
    shell:
        """
        Rscript {params.bin}/merge_patient_data.R \
            --tumor_pd1_dir {input.tumor_pd1_dir} \
            --risk_score {input.risk_score} \
            --immune_file {input.immune_file} \
            --output_file {output.out_file} \
            > {log} 2>&1
        """

## Plot the PC components vs PDL1 exp with risk scores from the univariate COX model ##

rule plot_PC_PDL1_expression:
    input:
        cox_summary_all = UNIVARIATE_COX_SUMMARY_ALL_FILTERED,
        tumor_main_dir = OUTPUT_DIR
    output:
        out_file = FIG_PC_PDL1_EXPRESSION
    log:
        "logs/plot_PC_PDL1_expression.log"
    message:
        "Generating a figure for the selected cancer types of PDL1 exp vs corresponding PC component"
    params:
        bin = config["bin"],
    shell:
        """
        Rscript {params.bin}/plot_PC_PDL1_exp.R \
            --cox_summary_all_cancers {input.cox_summary_all} \
            --tumor_dir {input.tumor_main_dir} \
            --output {output.out_file} \
            > {log} 2>&1
        """

## Plot the correlation of the PC components and the immune cell types ##

rule plot_PC_immune_correlations:
    input:
        cox_summary_all = UNIVARIATE_COX_SUMMARY_ALL_FILTERED,
        tumor_main_dir = OUTPUT_DIR
    output:
        out_file = FIG_PC_IMMUNE_CORRELATION
    log:
        "logs/plot_PC_immune_correlations.log"
    message:
        "Generating a figure of the correlation between the PC and immune cells for selected cancer types"
    params:
        bin = config["bin"],
    shell:
        """
        Rscript {params.bin}/plot_PC_immune_correlations.R \
            --cox_summary_all_cancers {input.cox_summary_all} \
            --tumor_dir {input.tumor_main_dir} \
            --output {output.out_file} \
            > {log} 2>&1
        """
## Perform the association analysis between the PD1 pathway heterogeneity scores and clinical features (categorical and numeric) ##

rule calculate_association_clinical_features_pd1_heterogeneity_scores:
    input:
        cox_res_file = UNIVARIATE_COX_SUMMARY_ALL_FILTERED,
        tumor_main_dir = OUTPUT_DIR
    output:
        results_pd1_groups = TUMOR_RESULTS_PD1_GROUPS,
        results_pd1_numeric = TUMOR_RESULTS_PD1_NUMERIC
    log:
        "logs/calculate_association_clinical_features_pd1_heterogeneity_scores.log"
    message:
        "Calculating PD1 heterogeneity clinical associations"
    params:
        bin = config["bin"],
    shell:
        """
        Rscript {params.bin}/clinical_associations_pd1.R \
            --cox_results_file {input.cox_res_file} \
            --tumor_main_dir {input.tumor_main_dir} \
            --results_pd1_groups {output.results_pd1_groups} \
            --results_pd1_numeric {output.results_pd1_numeric} \
            > {log} 2>&1
        """

## Plot the Associations between the PD1 pathway-based patient heterogeneity scores and clinical features

rule plot_association_clinical_features_pd1_heterogeneity_scores:
    input:
        cox_res_file = UNIVARIATE_COX_SUMMARY_ALL_FILTERED,
        results_pd1_groups = TUMOR_RESULTS_PD1_GROUPS,
        results_pd1_numeric = TUMOR_RESULTS_PD1_NUMERIC,
        cancer_color_file = CANCER_COLOR_FILE
    output:
        figure_pc_clin_associations = FIGURE_PC_CLIN_ASSOCIATIONS 
    log:
        "logs/plot_association_clinical_features_pd1_heterogeneity_scores.log"
    message:
        "Plot PD1 heterogeneity clinical associations"
    params:
        bin = config["bin"],
    shell:
        """
        Rscript {params.bin}/plot_clinical_associations_pd1.R \
            --cox_results_file {input.cox_res_file} \
            --results_pd1_groups {input.results_pd1_groups} \
            --results_pd1_numeric {input.results_pd1_numeric} \
            --cancer_color_file {input.cancer_color_file} \
            --pc_clinical_association_figure {output.figure_pc_clin_associations} \
            > {log} 2>&1

        """

## Plot the Associations between the PD1 pathway-based patient heterogeneity scores and individual clinical features

rule plot_individual_clinical_features_pd1_heterogeneity_scores:
    input:
        cox_res_file = UNIVARIATE_COX_SUMMARY_ALL_FILTERED,
        results_pd1_groups = TUMOR_RESULTS_PD1_GROUPS,
        results_pd1_numeric = TUMOR_RESULTS_PD1_NUMERIC,
        tumor_main_dir = OUTPUT_DIR
    output:
        figure_pc_individual_clin_associations = FIGURE_PC_INDIVIDUAL_CLIN_ASSOCIATIONS 
    log:
        "logs/plot_individual_clinical_features_pd1_heterogeneity_scores.log"
    message:
        "Plot PD1 heterogeneity clinical associations for individual features in the selected cancers"
    params:
        bin = config["bin"],
    shell:
        """
        Rscript {params.bin}/plot_individual_clinical_features_pd1_pathway.R \
            --cox_results_file {input.cox_res_file} \
            --results_pd1_groups {input.results_pd1_groups} \
            --results_pd1_numeric {input.results_pd1_numeric} \
            --tumor_main_dir {input.tumor_main_dir} \
            --output_figure_file {output.figure_pc_individual_clin_associations} \
            > {log} 2>&1

        """


## Run multivariate regularized cox regression on PDL1 edges ##  
rule run_regularized_cox:
    input:
        clin_file = OUTPUT_CANCER,
        tumor_pd1_dir = TUMOR_PD1_DIR
    output:
        out_file = OUTPUT_CANCER_COX
    message:
        "Running regularized cox for: {wildcards.cancer}"
    log:
        "logs/run_regularized_cox_{cancer}.log"
    params:
        bin = config["bin"],
        alpha = ALPHA,
        number_folds = NUMBER_FOLDS,
        number_cores = NUMBER_CORES,
        number_times = NUMBER_TIMES
    shell:
        """
        Rscript {params.bin}/cox_regression_tumor.R \
            --tumor_clin_file_path {input.clin_file} \
            --tumor_pd1_dir {input.tumor_pd1_dir} \
            --number_folds {params.number_folds} \
            --number_times {params.number_times} \
            --number_cores {params.number_cores} \
            --alpha {params.alpha} \
            --output {output.out_file} \
            > {log} 2>&1
        """
## Combine all multivarite cox results in one table ##  

rule combine_multivarite_cox_results:
    input:
        expand(OUTPUT_CANCER_COX, cancer=CANCER_TYPES)
    output:
        COX_RESULTS_ALL_MULTIVARIATE
    message:
        "Combining all multivariate Cox regression results into one table."
    log:
        "logs/combine_multivarite_cox_results.log"
    shell:
        """
        echo -e "cancer\\t$(head -n 1 {input[0]})" > {output}
        for file in {input}; do
            cancer=$(basename $(dirname $(dirname $file)))
            tail -n +2 $file | awk -v cancer=$cancer '{{print cancer"\\t"$0}}' >> {output}
        done
        echo "Finished combining results. Output written to: {output}" >> {log}
        """

## Create a circilar PDL1 plot with the selected TFs passing a threshold ##  

rule create_circular_pdl1_plot:
    input:
        cox_univariate_results_pd1_pathway = UNIVARIATE_COX_SUMMARY_ALL_FILTERED,
        cox_results_all_multivariate = COX_RESULTS_ALL_MULTIVARIATE,
        ppi_file = PPI_FILE,
        motif_file = MOTIF_FILE
    output:
        out_file = PDL1_CIRCULAR_PLOT
    log:
        "logs/create_circular_pdl1_plot_{threshold_cox}.log"
    message:
        "Making a PDL1 circular plot with threshold: {wildcards.threshold_cox}"
    params:
        bin = config["bin"],
    shell:
        """
        Rscript {params.bin}/circular_plot_PDL1.R \
            --cox_univariate_results_pd1_pathway {input.cox_univariate_results_pd1_pathway} \
            --cox_results_multivariate {input.cox_results_all_multivariate} \
            --ppi_file {input.ppi_file} \
            --motif_file {input.motif_file} \
            --threshold {wildcards.threshold_cox} \
            --output {output.out_file} \
            > {log} 2>&1
        """

## plot the comparison of the indegree and expression clusters for PRAD and UVM ##
rule plot_tsne_expression_indegree_and_uvm_prad_comparisons:
    input:
        datasets_to_plot_cola_clusters = DATASETS_TO_PLOT_COLA_CLUSTERS,
        tsne_file_expression = TSNE_DATA_EXPRESSION,
        tsne_file_indegree = TSNE_DATA_INDEGREE,
        cancer_color_file = CANCER_COLOR_FILE
    output:
        ouput_figure = FIGURE_TSNE_ALL_CANCERS_UVM_PRAD_CLUSTERS
    log:
        "logs/plot_tsne_expression_indegree_and_uvm_prad_comparisons.log"
    params:
        bin = config["bin"],
    shell:
        """
        Rscript {params.bin}/plot_TSNE_all_cancers.R \
            --datasets_to_plot_cola_clusters {input.datasets_to_plot_cola_clusters} \
            --tsne_data_expression {input.tsne_file_expression} \
            --tsne_data_indegree {input.tsne_file_indegree} \
            --cancer_color_file {input.cancer_color_file} \
            --ouput_figure_file {output.ouput_figure} \
            > {log} 2>&1

        """
## plot PD1 summary table ##
rule plot_PD1_summary_table:
    input:
        summary_table_PD1 = SUMMARY_TABLE_PD1,
    output:
        output_html_file = OUTPUT_HTML_TABLE_PD1
    log:
        "logs/plot_PD1_summary_table.log"
    params:
        bin = config["bin"],
    shell:
        """
        Rscript {params.bin}/plot_PD1_summary_table.R \
            --summary_table_PD1 {input.summary_table_PD1} \
            --output_html_file {output.output_html_file} \
            > {log} 2>&1
        """
