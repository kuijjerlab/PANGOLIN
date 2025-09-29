# cola consensus clustering on gene indegree and expression data ##
rule run_cola_clustering:
    input:
        expression_file = EXPRESSION_PANDA_FILE,
        samples_file = SAMPLES_WITH_CANCER_FILE,
        indegree_dir = INPUT_CANCER_INDEGREE_DIR
    output:
        best_k = BEST_K_COLA_OUTPUT_CANCER_DATATYPE,
        results = COLA_RESULTS_OUTPUT_CANCER_DATATYPE,
        collected_plots = COLA_COLLECTED_PLOTS_CANCER_DATATYPE,
        tsne_pdf = COLA_TSNE_OUTPUT_CANCER_DATATYPE,
        membership = COLA_MEMBERSHIP_OUTPUT_CANCER_DATATYPE,
        statistics = COLA_STATISTICS_OUTPUT_CANCER_DATATYPE,
        classes = COLA_CLASSES_OUTPUT_CANCER_DATATYPE,
        partition_pdf = COLA_PARTITION_OUTPUT_CANCER_DATATYPE,
       
    params:
        bin = config["bin"],
        number_cores_cola = NUMBER_CORES_COLA,
        partition_method = PARTITION_METHOD,
        top_value_method = TOP_VALUE_METHOD,
        max_k = MAX_K
    log:
        "logs/cola_clustering_{cancer}_{datatype}.log"
    shell:
        """
        Rscript workflow/bin/cola_clustering.R \
            --tumor {wildcards.cancer} \
            --exp_file {input.expression_file} \
            --samples_file {input.samples_file} \
            --indegree_dir {input.indegree_dir} \
            --datatype {wildcards.datatype} \
            --number_cores {params.number_cores_cola} \
            --top_value_method {params.top_value_method} \
            --partition_method {params.partition_method} \
            --max_k {params.max_k} \
            --output_best_k {output.best_k} \
            --output_results {output.results} \
            --output_membership {output.membership} \
            --output_statistics {output.statistics} \
            --output_classes {output.classes} \
            --output_collected_plots {output.collected_plots} \
            --output_tsne_pdf {output.tsne_pdf} \
            --output_partition_pdf {output.partition_pdf} \
            > {log} 2>&1
        """

# select the best k across all cancer types for indegree and expression data ##
rule select_best_k_cola:
    input:
        indegree_best_k_files = BEST_K_COLA_INDGEGREE_FILES,
        expression_best_k_files = BEST_K_COLA_EXPRESSION_FILES

    output:
        best_k_cola_ind_file = BEST_K_COLA_IND,
        best_k_cola_exp_file = BEST_K_COLA_EXP
    log:
        "logs/select_best_k_cola.log"
    shell:
        """
        Rscript workflow/bin/select_best_k_cola.R \
            --indegree_best_k_files "{input.indegree_best_k_files}" \
            --expression_best_k_files "{input.expression_best_k_files}" \
            --best_cola_k_indegree {output.best_k_cola_ind_file} \
            --best_cola_k_expression {output.best_k_cola_exp_file} \
            > {log} 2>&1
        """

rule plot_TSNE_cola_clusters:
    input:
        indegree_files = COLA_RESULTS_INDEGREE_FILES,
        expression_files = COLA_RESULTS_EXPRESSION_FILES,
        best_k_cola_ind_file = BEST_K_COLA_IND,
        best_k_cola_exp_file = BEST_K_COLA_EXP
    output:
        fig_tsne_indegree = FIG_TSNE_COLA_INDEGREE,
        fig_tsne_expression = FIG_TSNE_COLA_EXPRESSION
    log:
        "logs/plot_tsne_cola_clusters.log"
    params:
        bin = config["bin"],
    shell:
        """
        Rscript {params.bin}/TSNE_plot_cola_clusters.R \
            --indegree_files "{input.indegree_files}" \
            --expression_files "{input.expression_files}" \
            --best_cola_k_indegree {input.best_k_cola_ind_file} \
            --best_cola_k_expression {input.best_k_cola_exp_file} \
            --figure_TSNE_indegree {output.fig_tsne_indegree} \
            --figure_TSNE_expression {output.fig_tsne_expression} \
            > {log} 2>&1
        """
## Sanky plot comparing the indegree and expression clusters for each cancer type ##

rule plot_SANKEY_cola_clusters:
    input:
        indegree_files = COLA_RESULTS_INDEGREE_FILES,
        expression_files = COLA_RESULTS_EXPRESSION_FILES,
        best_k_cola_ind_file = BEST_K_COLA_IND,
        best_k_cola_exp_file = BEST_K_COLA_EXP
    output:
        fig_sankey_plot = FIG_SANKEY,
        selected_cola_ind_clusters_file = SELECTED_CLUSTERS_COLA_IND,
        selected_cola_exp_clusters_file = SELECTED_CLUSTERS_COLA_EXP,
        datasets_to_plot_cola_clusters = DATASETS_TO_PLOT_COLA_CLUSTERS
    log:
        "logs/plot_sankey_cola_clusters.log"
    message:
        "Plotting sankey plot comparing indegree and expression clusters"
    params:
        bin = config["bin"],
    shell:
        """
        Rscript {params.bin}/cola_clusters_sanky_plots.R \
            --indegree_files "{input.indegree_files}" \
            --expression_files "{input.expression_files}" \
            --best_cola_k_indegree {input.best_k_cola_ind_file} \
            --best_cola_k_expression {input.best_k_cola_exp_file} \
            --clusters_indegree {output.selected_cola_ind_clusters_file} \
            --clusters_expression {output.selected_cola_exp_clusters_file} \
            --datasets_to_plot_cola_clusters {output.datasets_to_plot_cola_clusters} \
            --figure_sanky {output.fig_sankey_plot} \
            > {log} 2>&1
        """
## Extract cola clusters for each cancer type and write to a separate files ##
rule save_final_cola_clusters_per_tumor:
    input:
        cluster_file_expression = SELECTED_CLUSTERS_COLA_EXP,
        cluster_file_indegree = SELECTED_CLUSTERS_COLA_IND,
    output:
        cluster_file_exp_per_cancer = OUTPUT_CLUSTERS_PER_TUMOR_EXP,
        cluster_file_ind_per_cancer = OUTPUT_CLUSTERS_PER_TUMOR_IND
    log:
        "logs/save_final_cola_clusters_per_tumor_{cancer}.log"
    params:
        bin = config["bin"]
    shell:
        """
        Rscript {params.bin}/extract_cola_clusters_per_tumor.R \
            --cluster_file_expression {input.cluster_file_expression} \
            --cluster_file_indegree {input.cluster_file_indegree} \
            --tumor {wildcards.cancer} \
            --cluster_expression_per_tumor {output.cluster_file_exp_per_cancer} \
            --cluster_indegree_per_tumor {output.cluster_file_ind_per_cancer} \
            > {log} 2>&1
        """

## Run univariate COX regression comparing cola clusters (expression and indegree) ##
rule run_univariate_cox_cola_clusters:
    input:
        clin_file = TUMOR_CLIN_FILE,
        cluster_file_exp_per_cancer = OUTPUT_CLUSTERS_PER_TUMOR_EXP,
        cluster_file_ind_per_cancer = OUTPUT_CLUSTERS_PER_TUMOR_IND
    output:
        out_file_summary = OUTPUT_CANCER_UNIVARIATE_COX_COLA_CLUSTERS
    log:
        "logs/run_univariate_cox_cola_clusters_{cancer}.log"
    message:
        "Running univariate Cox model comparing cola clusters for: {wildcards.cancer}"
    params:
        bin = config["bin"]
    shell:
        """
        Rscript {params.bin}/cola_clusters_survival.R \
            --clinical_file_tumor {input.clin_file} \
            --cluster_file_expression {input.cluster_file_exp_per_cancer} \
            --cluster_file_indegree {input.cluster_file_ind_per_cancer} \
            --output_file {output.out_file_summary} \
            > {log} 2>&1
        """
## Combine all COX regression results for comparing cola clusters (expression and indegree) ##
rule combine_univariate_cox_cola_clusters_results:
    input:
        expand(OUTPUT_CANCER_UNIVARIATE_COX_COLA_CLUSTERS, cancer=CANCER_TYPES)
    output:
        OUTPUT_CANCER_UNIVARIATE_COX_COLA_CLUSTERS_ALL
    log:
        "logs/combine_univariate_cox_cola_clusters_results.log"
    message:
        "Combining all univariate cox results (comparing cola clusters) into one table."
    shell:
        """
        echo -e "cancer\\t$(head -n 1 {input[0]})" > {output}
        for file in {input}; do
            cancer=$(basename $(dirname $(dirname $file)))
            tail -n +2 $file | awk -v cancer=$cancer '{{print cancer"\\t"$0}}' >> {output}
        done
        echo "Finished combining results. Output written to: {output}" >> {log}
        """

## Create a plot with COX results for cola clusters (expression and indegree) ##
rule plot_univariate_cox_cola_clusters_results:
    input:
        cox_cola_clusters_results = OUTPUT_CANCER_UNIVARIATE_COX_COLA_CLUSTERS_ALL,
        cancer_color_file = CANCER_COLOR_FILE
    output:
        fig_cox_cola_clusters = FIG_COX_COLA_CLUSTERS
    log:
        "logs/plot_univariate_cox_cola_clusters_results.log"
    message:
        "Plotting all univariate cox results (comparing cola clusters)"
    params:
        bin = config["bin"]
    shell:
        """
        Rscript {params.bin}/plot_cox_cola_clusters.R \
            --cox_results_cluster_file {input.cox_cola_clusters_results} \
            --cancer_color_file {input.cancer_color_file} \
            --output_file {output.fig_cox_cola_clusters} \
            > {log} 2>&1
        """
