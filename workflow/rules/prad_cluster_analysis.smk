## Create a survival plot for PRAD indegree cola clusters (k=4) ##
rule plot_PRAD_clusters_survival:
    input:
        prad_clin_file = PRAD_CLIN_FILE,
        prad_cluster_file_ind = CLUSTER_INDEGREE_PRAD,
        prad_indegree_file = PRAD_IND_FILE,
        expression_file = EXPRESSION_PANDA_FILE,
        samples_file = SAMPLES_WITH_CANCER_FILE,
        gmt_file = GMT_FILE
    output:
        fig_prad_survival = FIG_PRAD_SURVIVAL,
        fig_fgsea_prad = FIG_FGSEA_PRAD 
    log:
         "logs/plot_prad_clusters_survival.log"
    message:
        "Pltting the survival plot and fgsea results for PRAD indegree cola clusters"
    params:
        bin = config["bin"]
    shell:
        """
        Rscript {params.bin}/plot_prad_clusters_survival.R \
            --prad_clin_file_path {input.prad_clin_file} \
            --prad_cluster_file_indegree {input.prad_cluster_file_ind} \
            --indegree_file {input.prad_indegree_file} \
            --exp_file {input.expression_file} \
            --samples_file {input.samples_file} \
            --gmt_file {input.gmt_file} \
            --output_survival_plot {output.fig_prad_survival} \
            --output_fgsea_plot {output.fig_fgsea_prad} \
            > {log} 2>&1
        """
