## Extract clinical data for each cancer type ##
rule extract_clinical_data:
    input:
        clin_file = CLINICAL_FILE
    output:
        out_file = OUTPUT_CANCER
    log:
        "logs/extract_clinical_data_{cancer}.log"
    message:
        "Extracting clinical data for: {wildcards.cancer}"
    params:
        bin = config["bin"]
    shell:
        """
        Rscript {params.bin}/extract_clinical_data.R \
            --clinical {input.clin_file} \
            --tumor {wildcards.cancer} \
            --output {output.out_file} \
            > {log} 2>&1
        """

## Run T-SNE on expression and gene indegree data for all cancers ##
rule run_tsne:
    input:
        expression_file = EXPRESSION_PANDA_FILE,
        samples_file = SAMPLES_WITH_CANCER_FILE,
        tumor_main_dir = OUTPUT_DIR
    output:
        out_file_expression = TSNE_DATA_EXPRESSION,
        out_file_indegree = TSNE_DATA_INDEGREE
    log:
        "logs/run_tsne.log"
    message:
        "Running T-SNE on the indegree and expression"
    params:
        bin = config["bin"]
    shell:
        """
        Rscript {params.bin}/run_tsne.R \
            --exp_file {input.expression_file} \
            --samples_file {input.samples_file} \
            --tumor_dir {input.tumor_main_dir} \
            --output_file_expression {output.out_file_expression} \
            --output_file_indegree {output.out_file_indegree} \
            > {log} 2>&1
        """

## Create cancer legend for the cancer types ##
rule create_cancer_legend:
    input:
        cancer_color_file = CANCER_COLOR_FILE
    output:
        cancer_legend_pdf = CANCER_LEGEND_PDF
    log:
        "logs/create_cancer_legend.log"
    message:
        "Creating cancer legend"
    params:
        bin = config["bin"],
        fig_dir = FIG_DIR
    shell:
        """
        Rscript {params.bin}/create_cancer_legend.R \
            --cancer_color_file {input.cancer_color_file} \
            --output_file {output.cancer_legend_pdf} \
            > {log} 2>&1
        """
# filter PORCUPINE results
rule filter_porcupine_results:
    input:
        porcupine_file = PORCUPINE_FILE
    output:
        filtered_porcupine_file = FILTERED_PORCUPINE_FILE
    log:
        "logs/filter_porcupine_results_{cancer}.log"
    message:
        "Running filtering of porcupine results for: {wildcards.cancer}"
    params:
        bin = config["bin"]
    shell:
        """
        Rscript {params.bin}/preprocess_PORCUPINE_results.R \
            --porcupine_file_path {input.porcupine_file} \
            --filtered_porcupine_file_path {output.filtered_porcupine_file} \
            > {log} 2>&1
        """
# combine filtered PORCUPINE results
rule combine_porcupine_results:
    input:
        expand(FILTERED_PORCUPINE_FILE, cancer=CANCER_TYPES)
    output:
        PORCUPINE_RESULTS_ALL
    log:
        "logs/combine_porcupine_results.log"
    message:
        "Combining all filtered PORCUPINE results into one table."
    shell:
        """
       ( 
        echo -e "cancer\\t$(head -n 1 {input[0]})" > {output}
        for file in {input}; do
            cancer=$(basename $(dirname $(dirname $file)))
            tail -n +2 $file | awk -v cancer=$cancer '{{print cancer"\\t"$0}}' >> {output}
        done 
         )> {log} 2>&1
        """
# plot PORCUPINE results
rule plot_porcupine_results:
    input:
        porcupine_file_all = PORCUPINE_RESULTS_ALL,
        pathways_hierarchy_file = PATHWAYS_HIERARCHY_FILE,
        pathways_hsa_id_file = PATHWAYS_HSA_ID_FILE,
        list_of_pathways_file = LIST_PATHWAYS_FILE
    output:
        figure_pathway_intersection = FIG_PATHWAY_INTERSECTION,
        figure_shared_categories = FIG_SHARED_CATEGORIES
    log:
        "logs/plot_porcupine_results.log"
    message:
        "Creating plots for the PORCUPINE results" 
    params:
        bin = config["bin"]
    shell:
        """
        Rscript {params.bin}/plot_PORCUPINE_results.R \
            --pcp_results_all_cancers_file {input.porcupine_file_all} \
            --pathways_hierarchy_file {input.pathways_hierarchy_file} \
            --pathways_hsa_id_file {input.pathways_hsa_id_file} \
            --list_of_pathways_file {input.list_of_pathways_file} \
            --figure_pathway_intersection {output.figure_pathway_intersection} \
            --figure_shared_categories {output.figure_shared_categories} \
            > {log} 2>&1
        """
