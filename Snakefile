## PANGOLIN:
## A comprehensive Snakemake pipeline for TCGA analysis of gene regulatory networks ##

## How to run this pipeline:
## 1. Dry run (check workflow): snakemake --cores 1 -np
## 2. Load required modules: 
##    - module load snakemake/7.23.1-foss-2022a
##    - module load R-bundle-Bioconductor/3.15-foss-2022a-R-4.2.1
## 3. Source conda: source ~/.bashrc  # conda version 24.1.2
## 4. Clean conda cache: rm -rf .snakemake/conda
## 5. Execute pipeline: snakemake --use-conda --conda-frontend conda --cores 1

## PREREQUISITES:
## Before running, ensure you have:
## 1. PySNAIL package for normalization: 
##    git clone git@github.com:kuijjerlab/PySNAIL.git (place in envs/)
## 2. MBatch conda environment created from envs/mbatch.yaml:
##    conda env create -f envs/mbatch.yaml
## 3. All input data files specified in config.yaml

## WORKFLOW OVERVIEW:
## This pipeline performs comprehensive pan-cancer analysis including:
## - Batch effect detection and visualization using MBatch
## - Gene expression normalization with PySNAIL/qsmooth
## - Consensus clustering analysis with COLA
## - Survival analysis with Cox regression models  
## - Pathway enrichment analysis with PORCUPINE
## - t-SNE dimensionality reduction and visualization
## - Clinical associations and immune infiltration analysis




## Libraries
import os 
import sys
import glob
from pathlib import Path
import time

## Config

global CONFIG_PATH
CONFIG_PATH = "config.yaml"
configfile: CONFIG_PATH



## Directories ##
OUTPUT_DIR = config["output_dir"]
FIG_DIR = config["fig_dir"]

## From config ##
CANCER_TYPES = config["cancer_types"]
BATCH_FILE = config["batch_file"]
DATATYPES = config["datatypes"]
CLINICAL_FILE = config["clinical_file"]
CLINICAL_FILE_RDATA = config["clinical_file_rdata"]
CANCER_COLOR_FILE = config["cancer_color_file"]
PPI_FILE = config["ppi_file"]
MOTIF_FILE = config["motif_file"]
EXPRESSION_FILE = config["expression_file"]
SAMPLES_FILE = config["samples_file"]
IMMUNE_FILE = config["immune_file"]
GMT_FILE = config["gmt_file"]
PATHWAYS_HIERARCHY_FILE = config["pathways_hierarchy_file"]
PATHWAYS_HSA_ID_FILE = config["pathways_hsa_id_file"]
LIST_PATHWAYS_FILE = config["list_of_pathways_file"]

#####
TUMOR_CLIN_FILE = os.path.join(OUTPUT_DIR, "{cancer}", "clinical", "curated_clinical_{cancer}.txt")
PORCUPINE_FILE = os.path.join(OUTPUT_DIR, "{cancer}", "porcupine", "pcp_results_with_variance_{cancer}.txt")
TUMOR_PD1_DIR = os.path.join(OUTPUT_DIR, "{cancer}", "pd1_data")
TUMOR_PATHWAYS_MAPPING_PATH = os.path.join(OUTPUT_DIR, "{cancer}", "porcupine", "individual_scores_{cancer}.RData")
TSNE_DIR = os.path.join(OUTPUT_DIR, "tsne_results")
INPUT_CANCER_INDEGREE_DIR  = os.path.join(OUTPUT_DIR, "{cancer}", "indegrees_norm")

## Output directory for the downloaded GDC data ##
OUTPUT_GDC_DIR = os.path.join("data_all", "gdc_data")
OUTPUT_GDC_FILE = os.path.join(OUTPUT_GDC_DIR, "TCGA-{cancer}.RData")
MARKER_FILE =  os.path.join(OUTPUT_GDC_DIR, "{cancer}.SKIPPED.txt") # in case if we don't want to download the GDC data

# Output directory for the combined downloaded GDC data ##
EXPRESSION_DIR_GDC = os.path.join("data_all", "gdc_data_predownloaded") # directory for predownloaded GDC data
OUTPUT_DIR_DOWNLOAD_COMBINED = os.path.join("data_all", "combined_gdc_data")
OUTPUT_EXP_COMBINED_FILE = os.path.join(OUTPUT_DIR_DOWNLOAD_COMBINED, "hg38_STAR_counts.tsv")
GROUP_FILE = os.path.join(OUTPUT_DIR_DOWNLOAD_COMBINED, "hg38_sample_groups.tsv")
FEATURE_FILE = os.path.join(OUTPUT_DIR_DOWNLOAD_COMBINED, "hg38_features.RData")
PYSNAIL_NORMALIZED_FILE = os.path.join(OUTPUT_DIR_DOWNLOAD_COMBINED, "pysnail_normalized_STAR_counts.tsv")

# Output directory for the individual cancer expression files after normalization with PySNAIL ##

OUTPUT_DIR_PYSNAIL_CANCER = os.path.join("data_all", "pysnail_normalized_individual_cancer_expression")

# BATCH EFFECT ANALYSIS #
PYSNAIL_NORMALIZED_FILE_CANCER_SPECIFIC = os.path.join(OUTPUT_DIR_PYSNAIL_CANCER, "normalized_expression_TCGA-{cancer}.RData")
BATCH_DIR_ALL_CANCERS = os.path.join("data_all", "batch_analysis")
BATCH_DIR_CANCER = os.path.join(BATCH_DIR_ALL_CANCERS, "TCGA-{cancer}")
BATCH_EFFECT_PDF = os.path.join("figs", "MBatch_DSC.pdf")

# PANDA + LIONESS NETWORK INFERENCE #
NETWORKS_DIR = os.path.join("data_all", "networks")
# PANDA_NETWORK_FILE = os.path.join(NETWORKS_DIR, "panda_net.txt")


CANCER_LEGEND_PDF = os.path.join(FIG_DIR, "cancer_legend.pdf")
OUTPUT_CANCER = os.path.join(OUTPUT_DIR, "{cancer}", "clinical", "curated_clinical_{cancer}.txt")

## output files for filtered porcupine results ##
FILTERED_PORCUPINE_FILE = os.path.join(OUTPUT_DIR, "{cancer}", "porcupine", "pcp_results_filtered_{cancer}.txt")
PORCUPINE_RESULTS_ALL = os.path.join("data_all", "porcupine", "porcupine_results_all_cancers.txt")
FIG_PATHWAY_INTERSECTION = os.path.join(FIG_DIR, "pathways_intersection_pcp.pdf")
FIG_SHARED_CATEGORIES = os.path.join(FIG_DIR, "combined_figure_pcp_results_shared_categories.pdf")

## output directory for COLA-consesus clustering results for each cancer type ##
OUTPUT_CANCER_CONSENSUS_DIR = os.path.join(OUTPUT_DIR, "{cancer}", "consensus_clustering", "{datatype}")

## output directory for COLA-consesus clustering results for ALL cancer types ##
OUTPUT_ALL_CANCERS_CONSENSUS_DIR = os.path.join("data_all", "cola_consensus_clustering")
## output files  for COLA-consesus clustering results for ALL cancer types ##
BEST_K_COLA_EXP = os.path.join(OUTPUT_ALL_CANCERS_CONSENSUS_DIR, "best_k_cola_expression.txt")
BEST_K_COLA_IND = os.path.join(OUTPUT_ALL_CANCERS_CONSENSUS_DIR, "best_k_cola_indegree.txt")

## output files  for COLA-consesus clustering (individual to cluster assignment table) results for ALL cancer types ##
SELECTED_CLUSTERS_COLA_EXP = os.path.join(OUTPUT_ALL_CANCERS_CONSENSUS_DIR, "selected_clusters_expression.txt")
SELECTED_CLUSTERS_COLA_IND = os.path.join(OUTPUT_ALL_CANCERS_CONSENSUS_DIR, "selected_clusters_indegree.txt")
DATASETS_TO_PLOT_COLA_CLUSTERS = os.path.join(OUTPUT_ALL_CANCERS_CONSENSUS_DIR, "datasets_to_plot_cola_clusters.RData")

## output figures files  for COLA-consesus clustering results for ALL cancer types on TSNE ##
FIG_TSNE_COLA_INDEGREE = os.path.join(FIG_DIR, "TSNE_cola_clusters_indegree_all_cancers.pdf")
FIG_TSNE_COLA_EXPRESSION = os.path.join(FIG_DIR, "TSNE_cola_clusters_expression_all_cancers.pdf")
## output figures file for SANKEY plot comparing COLA clusters for indegree and expression ##
FIG_SANKEY = os.path.join(FIG_DIR, "sankey_plot_indegree_expression.pdf")

## output dir and files for saving cola clusters for each cancer type ##
OUTPUT_CLUSTERS_PER_TUMOR_IND = os.path.join(OUTPUT_DIR, "{cancer}", "final_clusters", "final_clusters_indegree_{cancer}.txt")
OUTPUT_CLUSTERS_PER_TUMOR_EXP = os.path.join(OUTPUT_DIR, "{cancer}", "final_clusters", "final_clusters_expression_{cancer}.txt")

## output file for cox univariate results comparing cola clustering results for each cancer type ##
OUTPUT_CANCER_UNIVARIATE_COX_COLA_CLUSTERS = os.path.join(OUTPUT_DIR, "{cancer}", "final_clusters", "cox_results_final_clusters_indegree_expression_{cancer}.txt")
OUTPUT_CANCER_UNIVARIATE_COX_COLA_CLUSTERS_ALL = os.path.join("data_all", "cox_results_all", "cox_results_final_clusters_indegree_expression_all.txt")

## output figure for cox univariate results comparing cola clustering results for all cancer types ##
FIG_COX_COLA_CLUSTERS = os.path.join(FIG_DIR, "cox_results_final_clusters_indegree_expression_all_cancers.pdf") 

## analysis of PRAD COLA CLUSTERS k = 4 ##
CLUSTER_INDEGREE_PRAD = os.path.join(OUTPUT_DIR, "PRAD", "final_clusters", "final_clusters_indegree_PRAD.txt")
PRAD_CLIN_FILE = os.path.join(OUTPUT_DIR, "PRAD", "clinical", "curated_clinical_PRAD.txt")
PRAD_PD1_DIR = os.path.join(OUTPUT_DIR, "PRAD", "pd1_data")
PRAD_IND_FILE = os.path.join(OUTPUT_DIR, "PRAD", "indegrees_norm", "indegree_norm_PRAD.RData")
FIG_PRAD_SURVIVAL = os.path.join(FIG_DIR, "PRAD_clusters_survival.pdf")
FIG_FGSEA_PRAD = os.path.join(FIG_DIR, "PRAD_clusters_fgsea.pdf")

## output directory for plot TSNE cola clusters ##
TSNE_DATA_DIR = os.path.join("data_all", "tsne_results")
TSNE_DATA_EXPRESSION = os.path.join(TSNE_DATA_DIR, "tsne_expression_all_cancers.txt")
TSNE_DATA_INDEGREE = os.path.join(TSNE_DATA_DIR, "tsne_expression_all_indegree.txt")

## Output Files for Univariate Cox on PD1-pathway based heterogeneity scores ##
OUTPUT_CANCER_UNIVARIATE_COX_SUMMARY = os.path.join(OUTPUT_DIR, "{cancer}", "cox", "{cancer}_PD1_pathway_cox_univariate_model_summary.txt")
OUTPUT_CANCER_UNIVARIATE_COX_PREDICTED_SCORES = os.path.join(OUTPUT_DIR, "{cancer}", "cox", "{cancer}_PD1_pathway_cox_univariate_predited_risk_scores.txt")
UNIVARIATE_COX_SUMMARY_ALL = os.path.join("data_all", "cox_results_all", "PD1_pathway_cox_univariate_model_summary_all.txt")
UNIVARIATE_COX_SUMMARY_ALL_FILTERED = os.path.join("data_all", "cox_results_all", "PD1_pathway_cox_univariate_model_summary_filtered.txt")
UNIVARIATE_COX_PREDICTED_SCORES_ALL = os.path.join("data_all", "cox_results_all", "PD1_pathway_cox_univariate_predited_risk_scores_all.txt")
FIG_PC_PDL1_EXPRESSION = os.path.join(FIG_DIR, "PDL1_exp_PC_component_HR.pdf")
FIG_PC_IMMUNE_CORRELATION = os.path.join(FIG_DIR, "PC_immune_correlations_cibersort.png")

TUMOR_RESULTS_PD1_GROUPS = os.path.join("data_all", "clinical_associations_PD1", "pd1_pathway_categorical_results.txt")
TUMOR_RESULTS_PD1_NUMERIC = os.path.join("data_all", "clinical_associations_PD1", "pd1_pathway_numeric_results.txt")
FIGURE_PC_CLIN_ASSOCIATIONS = os.path.join("figs", "PC_all_features_clin_associations.pdf")
FIGURE_PC_INDIVIDUAL_CLIN_ASSOCIATIONS = os.path.join("figs", "PC_individual_features_clin_associations.pdf")

## Output Files for multivariate regularized Cox on PDL1-edges ##
OUTPUT_CANCER_PD1_MAPPINGS  = os.path.join(OUTPUT_DIR, "{cancer}", "pd1_data", "pd1_individual_scores_norm_{cancer}.RData")
OUTPUT_CANCER_COX = os.path.join(OUTPUT_DIR, "{cancer}", "cox", "{cancer}_PDL1_cox_multivariate_res.txt")
COX_RESULTS_ALL_MULTIVARIATE = os.path.join("data_all", "cox_results_all", "PDL1_cox_multivarite_res_all.txt")
PDL1_CIRCULAR_PLOT = os.path.join(FIG_DIR, "circular_pdl1_plot_{threshold_cox}.pdf")
OUTPUT_PDL1_EXP_CANCER = os.path.join(OUTPUT_DIR, "{cancer}", "pd1_data", "pdl1_expression_{cancer}.txt")
OUTPUT_COMBINED_PATIENT_DATA_CANCER = os.path.join(OUTPUT_DIR, "{cancer}", "pd1_data", "combined_patient_data_{cancer}.txt")


# figure for TNSE plot for all cancers and also comparisons of cola clusters for indegree and expression for PRAD and UVM
FIGURE_TSNE_ALL_CANCERS_UVM_PRAD_CLUSTERS = os.path.join(FIG_DIR, "TSNE_all_cancers_indegree_expression_UVM_PRAD_clusters.pdf")

# producing summary table figure for PD1 pathway 
#input (manually created)
SUMMARY_TABLE_PD1 = os.path.join("data_all", "clinical_associations_PD1", "summary_table_PD1.txt")
#output
OUTPUT_HTML_TABLE_PD1 = os.path.join("figs", "summary_table_PD1.html")

## Parameters ##
ALPHA = config["alpha"]
NUMBER_FOLDS = config["number_folds"]
NUMBER_CORES = config["number_cores"]
NUMBER_TIMES = config["number_times"]
THESHOLD_COX = config["threshold_cox"]
GENE_ID = config["gene_id"]
NUMBER_CORES_COLA = config["number_cores_cola"]
PARTITION_METHOD = config["partition_method"]
TOP_VALUE_METHOD = config["top_value_method"]
MAX_K = config["max_k"]

# Rules ##
rule all:
    input:
        # expand(OUTPUT_GDC_FILE, cancer = CANCER_TYPES),
        # expand(MARKER_FILE, cancer = CANCER_TYPES),
        # OUTPUT_EXP_COMBINED_FILE, 
        # GROUP_FILE,
        # FEATURE_FILE,
        # PYSNAIL_NORMALIZED_FILE,
        # OUTPUT_DIR_PYSNAIL_CANCER,
        # expand(BATCH_DIR_CANCER, cancer = CANCER_TYPES),
        # BATCH_EFFECT_PDF
        NETWORKS_DIR
        # expand(OUTPUT_CANCER, cancer = CANCER_TYPES),
        # TSNE_DATA_EXPRESSION,
        # TSNE_DATA_INDEGREE,
        # CANCER_LEGEND_PDF,
        # expand(FILTERED_PORCUPINE_FILE, cancer = CANCER_TYPES),
        # PORCUPINE_RESULTS_ALL,
        # FIG_PATHWAY_INTERSECTION,
        # FIG_SHARED_CATEGORIES,
        # expand(OUTPUT_CANCER_CONSENSUS_DIR, cancer = CANCER_TYPES, datatype = DATATYPES),
        # BEST_K_COLA_EXP,
        # BEST_K_COLA_IND,
        # FIG_TSNE_COLA_INDEGREE,
        # FIG_TSNE_COLA_EXPRESSION,
        # FIG_SANKEY,
        # SELECTED_CLUSTERS_COLA_EXP,
        # SELECTED_CLUSTERS_COLA_IND,
        # expand(OUTPUT_CLUSTERS_PER_TUMOR_IND, cancer = CANCER_TYPES),
        # expand(OUTPUT_CLUSTERS_PER_TUMOR_EXP, cancer = CANCER_TYPES),
        # expand(OUTPUT_CANCER_UNIVARIATE_COX_COLA_CLUSTERS, cancer = CANCER_TYPES),
        # OUTPUT_CANCER_UNIVARIATE_COX_COLA_CLUSTERS_ALL,
        # FIG_COX_COLA_CLUSTERS,
        # FIG_PRAD_SURVIVAL,
        # FIG_FGSEA_PRAD, 
        # expand(OUTPUT_PDL1_EXP_CANCER, cancer = CANCER_TYPES),
        # expand(OUTPUT_CANCER_PD1_MAPPINGS, cancer = CANCER_TYPES),
        # expand(OUTPUT_CANCER_UNIVARIATE_COX_SUMMARY, cancer = CANCER_TYPES),
        # expand(OUTPUT_CANCER_UNIVARIATE_COX_PREDICTED_SCORES, cancer = CANCER_TYPES),
        # UNIVARIATE_COX_SUMMARY_ALL,
        # UNIVARIATE_COX_PREDICTED_SCORES_ALL,
        # UNIVARIATE_COX_SUMMARY_ALL_FILTERED,
        # expand(OUTPUT_COMBINED_PATIENT_DATA_CANCER, cancer = CANCER_TYPES),
        # FIG_PC_PDL1_EXPRESSION,
        # FIG_PC_IMMUNE_CORRELATION,
        # TUMOR_RESULTS_PD1_GROUPS,
        # TUMOR_RESULTS_PD1_NUMERIC,
        # FIGURE_PC_CLIN_ASSOCIATIONS,
        # FIGURE_PC_INDIVIDUAL_CLIN_ASSOCIATIONS,
        # expand(OUTPUT_CANCER_COX, cancer = CANCER_TYPES),
        # COX_RESULTS_ALL_MULTIVARIATE, 
        # expand(PDL1_CIRCULAR_PLOT, threshold_cox = THESHOLD_COX),
        # FIGURE_TSNE_ALL_CANCERS_UVM_PRAD_CLUSTERS,
        # OUTPUT_HTML_TABLE_PD1   

## Download GDC data for each cancer type ##    
## if the download_files is set to "YES" in the config file, then it will download the GDC data
## otherwise it will skip the download and create a marker file 
## and an empty output file indicating that the download was skipped
## The predownloaded GDC data is in the data_all/gdc_data_predownloaded directory

# rule download_gdc_data:
#     output:
#         out_file = OUTPUT_GDC_FILE,
#         marker = MARKER_FILE
#     message:
#         "Downloading GDC data for: {wildcards.cancer}"
#     params:
#         bin = config["bin"]
#     run:
#         if config["download_files"] == "YES":
#             shell(
#                 """
#                 Rscript {params.bin}/download_TCGA_expression.R \
#                     --tumor {wildcards.cancer} \
#                     --output_file {output.out_file}
#                 """
#             )
#             with open(output.marker, "w") as f:
#                 f.write("Downloaded\n")
#         else:
#             print(f"Skipping download for {wildcards.cancer}")
#             os.makedirs(os.path.dirname(output.marker), exist_ok=True)
#             with open(output.marker, "w") as f:
#                 f.write("Download skipped.\n")

# # Combine the data for all cancer types, including expression, groups and features ##
# rule combine_all_expression_data:
#     input:
#         exp_dir = EXPRESSION_DIR_GDC
#     output:
#         combined_expression_file = OUTPUT_EXP_COMBINED_FILE,
#         group_file = GROUP_FILE,
#         feature_file = FEATURE_FILE
#     message:
#         "Combining all expression data for all cancers"
#     params:
#         bin = config["bin"]
#     shell:
#         """
#         Rscript {params.bin}/combine_allTCGA_expression.R \
#             --expression_dir {input.exp_dir} \
#             --combined_expression_file {output.combined_expression_file} \
#             --group_file {output.group_file} \
#             --feature_file {output.feature_file}
#         """
# # SNAIL qsmooth normalization of the combined expression data #
# # to run this rule,you need to download the pysnail package in envs 
# # to download it git clone git@github.com:kuijjerlab/PySNAIL.git
# rule qsmooth_normalization:
#     input:
#         xprs = OUTPUT_EXP_COMBINED_FILE,
#         groups = GROUP_FILE
#     output:
#         norm = PYSNAIL_NORMALIZED_FILE
#     params:
#         threshold = 0.2,
#         bin = config["bin"]
#     conda:
#         "envs/pysnail.yaml"
#     shell:
#         """
#             python {params.bin}/normalize_with_pysnail.py {input.xprs} {input.groups} {output.norm} --threshold {params.threshold}
#         """

# PANDA/LIONESS network inference using netZooPy
rule run_panda_lioness:
    """
    Run PANDA and LIONESS network inference using netZooPy.
    PANDA builds an initial regulatory network and LIONESS estimates 
    sample-specific networks for each individual sample.
    """
    input:
        exp_file = EXPRESSION_FILE,
        motif_file = MOTIF_FILE,
        ppi_file = PPI_FILE
    output:
        network_dir = directory(NETWORKS_DIR)
    params:
        bin = config["bin"],
        start_sample = 1,
        end_sample = 10,  # adjust based on your sample count
        computing = "cpu",  # change to "gpu" if you have GPU support
        random_seed = 10,
        ncores = 10
    conda:
        "envs/netzoopy-local.yaml"
    shell:
        """
        set +u
        unset PYTHONPATH
        unset PYTHONHOME
        export PYTHONNOUSERSITE=1
        
        # Create networks directory
        mkdir -p {output.network_dir}
        
        python {params.bin}/run_panda_lioness.py \
            --exp_file {input.exp_file} \
            --motif_file {input.motif_file} \
            --ppi_file {input.ppi_file} \
            --output_dir {output.network_dir} \
            --start_sample {params.start_sample} \
            --end_sample {params.end_sample} \
            --computing {params.computing} \
            --random_seed {params.random_seed} \
            --ncores {params.ncores}
        """

# rule split_expression_by_cancer:
#     """
#     Split the large normalized expression file into cancer-specific files.
#     """
#     input:
#         expression_file = PYSNAIL_NORMALIZED_FILE,
#         group_file = GROUP_FILE
#     output:
#         output_directory = directory(OUTPUT_DIR_PYSNAIL_CANCER)
#     message:
#         "Splitting expression data into cancer-specific files and saving them"
#     params:
#         bin = config["bin"]
#     shell:
#         """
#         Rscript {params.bin}/save_normalized_exp_per_cancer.R \
#             --expression_file {input.expression_file} \
#             --group_file {input.group_file} \
#             --output_dir {output.output_directory}
#         """   
## Check the batch effect in the expression data ##
## This takes a while to run ##
rule analyze_batch_effect:
    input:
        expression_file = PYSNAIL_NORMALIZED_FILE_CANCER_SPECIFIC,
        group_file = GROUP_FILE,
        batch_file = BATCH_FILE,
        clinical_file = CLINICAL_FILE_RDATA 
    output:
        output_dir = directory(BATCH_DIR_CANCER)
    message:
        "Analyzing batch effect for: {wildcards.cancer}"   

    params:
        bin = config["bin"],
        nthreads = config["number_cores_mbatch"]
    conda:
        "mbatch_minimal"
    shell:
        """
        Rscript {params.bin}/batch_effect_analysis.R \
            --tumor_type {wildcards.cancer} \
            --expression_file {input.expression_file} \
            --group_file {input.group_file} \
            --batch_file {input.batch_file} \
            --clin_file {input.clinical_file} \
            --nthreads {params.nthreads} \
            --output_directory {output.output_dir}
        """

## Create a PDF figure for the batch effect analysis ##
rule plot_mbatch_results:
    input:
        batch_dir = BATCH_DIR_ALL_CANCERS
    output:
        pdf_file = BATCH_EFFECT_PDF
    message:
        "Creating batch effect figure"
    params:
        bin = config["bin"]
    shell:
        """
        Rscript {params.bin}/plot_mbatch_results.R \
            --batch_results_dir {input.batch_dir} \
            --output_file {output.pdf_file}
        """

## Extract clinical data for each cancer type ##
rule extract_clinical_data:
    input:
        clin_file = CLINICAL_FILE
    output:
        out_file = OUTPUT_CANCER
    message:
        "Extracting clinical data for: {wildcards.cancer}"
    params:
        bin = config["bin"]
    shell:
        """
        Rscript {params.bin}/extract_clinical_data.R \
            --clinical {input.clin_file} \
            --tumor {wildcards.cancer} \
            --output {output.out_file}
        """


## Run T-SNE on expression and gene indegree data for all cancers ##
rule run_tsne:
    input:
        expression_file = EXPRESSION_FILE,
        samples_file = SAMPLES_FILE,
        tumor_main_dir = OUTPUT_DIR
    output:
        out_file_expression = TSNE_DATA_EXPRESSION,
        out_file_indegree = TSNE_DATA_INDEGREE
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
            --output_file_indegree {output.out_file_indegree} 
            
        """

## Create cancer legend for the cancer types ##
rule create_cancer_legend:
    input:
        cancer_color_file = CANCER_COLOR_FILE
    output:
        CANCER_LEGEND_PDF
    params:
        bin = config["bin"],
        fig_dir = FIG_DIR
    shell:
        """
        Rscript {params.bin}/create_cancer_legend.R \
            --cancer_color_file {input.cancer_color_file} \
            --figure_dir {params.fig_dir}
        """
# filter PORCUPINE results
rule filter_porcupine_results:
    input:
        porcupine_file = PORCUPINE_FILE
    output:
        filtered_porcupine_file = FILTERED_PORCUPINE_FILE
    message:
        "Running filtering of porcupine results for: {wildcards.cancer}" 
    params:
        bin = config["bin"]
    shell:
        """
        Rscript {params.bin}/preprocess_PORCUPINE_results.R \
            --porcpupine_file_path {input.porcupine_file} \
            --filtered_porcupine_file_path {output.filtered_porcupine_file}
        """
# combine filtered PORCUPINE results
rule combine_porcupine_results:
    input:
        expand(FILTERED_PORCUPINE_FILE, cancer=CANCER_TYPES)
    output:
        PORCUPINE_RESULTS_ALL
    message:
        "Combining all filtered PORCUPINE results into one table."
    shell:
        """
        echo -e "cancer\\t$(head -n 1 {input[0]})" > {output}
        for file in {input}; do
            cancer=$(basename $(dirname $(dirname $file)))
            tail -n +2 $file | awk -v cancer=$cancer '{{print cancer"\\t"$0}}' >> {output}
        done
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
            --figure_shared_categories {output.figure_shared_categories}
        """

# cola consensus clustering on gene indegree and expression data ##
rule run_cola_clustering:
    input:
        expression_file = EXPRESSION_FILE,
        samples_file = SAMPLES_FILE,
        indegree_dir = INPUT_CANCER_INDEGREE_DIR
    output:
        output_cancer_consensus_dir = directory(OUTPUT_CANCER_CONSENSUS_DIR)
    message:
        "Running cola clustering on: {wildcards.cancer} with datatype: {wildcards.datatype}"
    params:
        bin = config["bin"],
        number_cores_cola = NUMBER_CORES_COLA,
        partition_method = PARTITION_METHOD,
        top_value_method = TOP_VALUE_METHOD,
        max_k = MAX_K
    shell:
        """
        Rscript {params.bin}/cola_clustering.R \
            --tumor {wildcards.cancer} \
            --exp_file {input.expression_file} \
            --samples_file {input.samples_file} \
            --indegree_dir {input.indegree_dir} \
            --datatype {wildcards.datatype} \
            --number_cores {params.number_cores_cola} \
            --top_value_method {params.top_value_method} \
            --partition_method {params.partition_method} \
            --max_k {params.max_k} \
            --output {output.output_cancer_consensus_dir}
        """

# select best K suggested by cola for each cancer type (indegree and expression) ##
rule select_best_k_cola:
    input:
        tumor_main_dir = OUTPUT_DIR
    output:
        best_k_cola_ind_file = BEST_K_COLA_IND,
        best_k_cola_exp_file = BEST_K_COLA_EXP
    message:
        "Selecting best number of clusters for each cancer"
    params:
        bin = config["bin"],
    shell:
        """
        Rscript {params.bin}/select_best_k_cola.R \
            --tumor_dir {input.tumor_main_dir} \
            --best_cola_k_indegree {output.best_k_cola_ind_file} \
            --best_cola_k_expression {output.best_k_cola_exp_file} \
        """
# plot TSNE with the selected cola clusters for each cancer type (indegree and expression) ##
rule plot_TSNE_cola_clusters:
    input:
        tumor_main_dir = OUTPUT_DIR,
        best_k_cola_ind_file = BEST_K_COLA_IND,
        best_k_cola_exp_file = BEST_K_COLA_EXP,
    output:
        fig_tsne_indegree = FIG_TSNE_COLA_INDEGREE,
        fig_tsne_expression = FIG_TSNE_COLA_EXPRESSION
    message:
        "Plotting T-SNE with cola clusters for all cancer types"
    params:
        bin = config["bin"],
    shell:
        """
        Rscript {params.bin}/TSNE_plot_cola_clusters.R \
            --tumor_dir {input.tumor_main_dir} \
            --best_cola_k_indegree {input.best_k_cola_ind_file} \
            --best_cola_k_expression {input.best_k_cola_exp_file} \
            --figure_TSNE_indegree {output.fig_tsne_indegree} \
            --figure_TSNE_expression {output.fig_tsne_expression} 
        """


## Sanky plot comparing the indegree and expression clusters for each cancer type ##
rule plot_SANKEY_cola_clusters:
    input:
        tumor_main_dir = OUTPUT_DIR,
        best_k_cola_ind_file = BEST_K_COLA_IND,
        best_k_cola_exp_file = BEST_K_COLA_EXP
    output:
        fig_sankey_plot = FIG_SANKEY,
        selected_cola_ind_clusters_file = SELECTED_CLUSTERS_COLA_IND,
        selected_cola_exp_clusters_file = SELECTED_CLUSTERS_COLA_EXP,
        datasets_to_plot_cola_clusters = DATASETS_TO_PLOT_COLA_CLUSTERS
    message:
        "Plotting sankey plot comparing indegree and expression clusters"
    params:
        bin = config["bin"],
    shell:
        """
        Rscript {params.bin}/cola_clusters_sanky_plots.R \
            --tumor_dir {input.tumor_main_dir} \
            --best_cola_k_indegree {input.best_k_cola_ind_file} \
            --best_cola_k_expression {input.best_k_cola_exp_file} \
            --clusters_indegree {output.selected_cola_ind_clusters_file} \
            --clusters_expression {output.selected_cola_exp_clusters_file} \
            --datasets_to_plot_cola_clusters {output.datasets_to_plot_cola_clusters} \
            --figure_sanky {output.fig_sankey_plot}  
        """

## Extract cola clusters for each cancer type and write to a separate files ##
rule save_final_cola_clusters_per_tumor:
    input:
        cluster_file_expression = SELECTED_CLUSTERS_COLA_EXP,
        cluster_file_indegree = SELECTED_CLUSTERS_COLA_IND,
    output:
        cluster_file_exp_per_cancer = OUTPUT_CLUSTERS_PER_TUMOR_EXP,
        cluster_file_ind_per_cancer = OUTPUT_CLUSTERS_PER_TUMOR_IND
    params:
        bin = config["bin"]
    shell:
        """
        Rscript {params.bin}/extract_cola_clusters_per_tumor.R \
            --cluster_file_expression {input.cluster_file_expression} \
            --cluster_file_indegree {input.cluster_file_indegree} \
            --tumor {wildcards.cancer} \
            --cluster_expression_per_tumor {output.cluster_file_exp_per_cancer} \
            --cluster_indegree_per_tumor {output.cluster_file_ind_per_cancer}
        """

## Run univariate COX regression comparing cola clusters (expression and indegree) ##
rule run_univariate_cox_cola_clusters:
    input:
        clin_file = TUMOR_CLIN_FILE,
        tumor_pd1_dir = TUMOR_PD1_DIR,
        cluster_file_exp_per_cancer = OUTPUT_CLUSTERS_PER_TUMOR_EXP,
        cluster_file_ind_per_cancer = OUTPUT_CLUSTERS_PER_TUMOR_IND
    output:
        out_file_summary = OUTPUT_CANCER_UNIVARIATE_COX_COLA_CLUSTERS
    message:
        "Running univariate Cox model comparing cola clusters for: {wildcards.cancer}"
    params:
        bin = config["bin"]
    shell:
        """
        Rscript {params.bin}/cola_clusters_survival.R \
            --clinical_file_tumor {input.clin_file} \
            --tumor_pd1_directory {input.tumor_pd1_dir} \
            --cluster_file_expression {input.cluster_file_exp_per_cancer} \
            --cluster_file_indegree {input.cluster_file_ind_per_cancer} \
            --output_file {output.out_file_summary}
        """
## Combine all COX regression results for comparing cola clusters (expression and indegree) ##
rule combine_univariate_cox_cola_clusters_results:
    input:
        expand(OUTPUT_CANCER_UNIVARIATE_COX_COLA_CLUSTERS, cancer=CANCER_TYPES)
    output:
        OUTPUT_CANCER_UNIVARIATE_COX_COLA_CLUSTERS_ALL
    message:
        "Combining all univariate cox results (comparing cola clusters) into one table."
    shell:
        """
        echo -e "cancer\\t$(head -n 1 {input[0]})" > {output}
        for file in {input}; do
            cancer=$(basename $(dirname $(dirname $file)))
            tail -n +2 $file | awk -v cancer=$cancer '{{print cancer"\\t"$0}}' >> {output}
        done
        """

## Create a plot with COX results for cola clusters (expression and indegree) ##
rule plot_univariate_cox_cola_clusters_results:
    input:
        cox_cola_clusters_results = OUTPUT_CANCER_UNIVARIATE_COX_COLA_CLUSTERS_ALL,
        cancer_color_file = CANCER_COLOR_FILE
    output:
        fig_cox_cola_clusters = FIG_COX_COLA_CLUSTERS
    message:
        "Pltting all univariate cox results (comparing cola clusters)"
    params:
        bin = config["bin"]
    shell:
        """
        Rscript {params.bin}/plot_cox_cola_clusters.R \
            --cox_results_cluster_file {input.cox_cola_clusters_results} \
            --cancer_color_file {input.cancer_color_file} \
            --output_file {output.fig_cox_cola_clusters} \
        """


## Create a survival plot for PRAD indegree cola clusters (k=4) ##
rule plot_PRAD_clusters_survival:
    input:
        prad_clin_file = PRAD_CLIN_FILE,
        prad_pd1_dir = PRAD_PD1_DIR,
        prad_cluster_file_ind = CLUSTER_INDEGREE_PRAD,
        prad_indegree_file = PRAD_IND_FILE,
        expression_file = EXPRESSION_FILE,
        samples_file = SAMPLES_FILE, 
        gmt_file = GMT_FILE
    output:
        fig_prad_survival = FIG_PRAD_SURVIVAL,
        fig_fgsea_prad = FIG_FGSEA_PRAD 
    message:
        "Pltting the survival plot and fgsea results for PRAD indegree cola clusters"
    params:
        bin = config["bin"]
    shell:
        """
        Rscript {params.bin}/plot_prad_clusters_survival.R \
            --prad_clin_file_path {input.prad_clin_file} \
            --prad_pd1_dir {input.prad_pd1_dir} \
            --prad_cluster_file_indegree {input.prad_cluster_file_ind} \
            --indegree_file {input.prad_indegree_file} \
            --exp_file {input.expression_file} \
            --samples_file {input.samples_file} \
            --gmt_file {input.gmt_file} \
            --output_survival_plot {output.fig_prad_survival} \
            --output_fgsea_plot {output.fig_fgsea_prad}
        """


## Extract PDL1 gene expression for each cancer type ##

rule extract_PDL1_gene_expression:
    input:
        expression_file = EXPRESSION_FILE,
        samples_file = SAMPLES_FILE
    output:
        out_file = OUTPUT_PDL1_EXP_CANCER
    message:
        "Extracting expression data for PDL1: {wildcards.cancer}"
    params:
        bin = config["bin"],
        gene_id = GENE_ID # this is how it shouldbe

    shell:
        """
        Rscript {params.bin}/extract_PDL1_expression.R \
            --exp_file {input.expression_file} \
            --samples_file {input.samples_file} \
            --tumor {wildcards.cancer} \
            --gene_id {params.gene_id} \
            --output {output.out_file}
        """

## Extract PD1 pathway-based individual mappings (heterogeneity scores) ##        
rule extract_pd1_pathway_individual_mappings:
    input:
       tumor_pathways_mapping_path = TUMOR_PATHWAYS_MAPPING_PATH
    output:
        out_file = OUTPUT_CANCER_PD1_MAPPINGS
    params:
        bin = config["bin"]
    message:
        "Extracting PD1 pathway individual mappings: {wildcards.cancer}"
    shell:
        """
        Rscript {params.bin}/extract_pd1_pathway_individual_scores.R \
            --tumor_pathways_mapping_path {input.tumor_pathways_mapping_path} \
            --output {output.out_file}
        """


## Run univatiate cox model on pd1-pathway based scores for each cancer ##        
rule run_univariate_cox_pd1_pathway:
    input:
        clin_file = TUMOR_CLIN_FILE,
        tumor_pd1_dir = TUMOR_PD1_DIR
    output:
        out_file_summary = OUTPUT_CANCER_UNIVARIATE_COX_SUMMARY,
        out_file_predicted_scores = OUTPUT_CANCER_UNIVARIATE_COX_PREDICTED_SCORES
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
            --cox_predicted_risk {output.out_file_predicted_scores}
        """
rule combine_pd1_pathway_univariate_summary_results:
    input:
        expand(OUTPUT_CANCER_UNIVARIATE_COX_SUMMARY, cancer=CANCER_TYPES)
    output:
        UNIVARIATE_COX_SUMMARY_ALL
    message:
        "Combining all univariate Cox regression results into one table."
    shell:
        """
        echo -e "cancer\\t$(head -n 1 {input[0]})" > {output}
        for file in {input}; do
            cancer=$(basename $(dirname $(dirname $file)))
            tail -n +2 $file | awk -v cancer=$cancer '{{print cancer"\\t"$0}}' >> {output}
        done
        """

rule combine_pd1_pathway_univariate_prediction_scores:
    input:
        expand(OUTPUT_CANCER_UNIVARIATE_COX_PREDICTED_SCORES, cancer=CANCER_TYPES)
    output:
        UNIVARIATE_COX_PREDICTED_SCORES_ALL
    message:
        "Combining all univariate Cox prediction risk results into one table."
    shell:
        """
        echo -e "cancer\\t$(head -n 1 {input[0]})" > {output}
        for file in {input}; do
            cancer=$(basename $(dirname $(dirname $file)))
            tail -n +2 $file | awk -v cancer=$cancer '{{print cancer"\\t"$0}}' >> {output}
        done
        """

## Clean the univatiate cox model on pd1-pathway results based of the PORCUPINE results ##        
rule clean_univariate_cox_pd1_pathway:
    input:
        univarite_cox_summary_all = UNIVARIATE_COX_SUMMARY_ALL,
        porcupine_results_all_filtered = PORCUPINE_RESULTS_ALL
    output:
        univarite_cox_summary_all_filtered = UNIVARIATE_COX_SUMMARY_ALL_FILTERED
    params:
        bin = config["bin"]
    message:
        "Cleaning the univariate Cox model on pd1-pathway results based on the PORCUPINE results"
    shell:
        """
        Rscript {params.bin}/clean_pd1_pathway_cox_results_by_porcupine.R \
            --cox_summary_all_cancers {input.univarite_cox_summary_all} \
            --porcupine_filtered_results {input.porcupine_results_all_filtered} \
            --cox_summary_all_cancers_filtered {output.univarite_cox_summary_all_filtered} 
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
            --output_file {output.out_file}
        """

## Plot the PC components vs PDL1 exp with risk scores from the univariate COX model ##

rule plot_PC_PDL1_expression:
    input:
        cox_summary_all = UNIVARIATE_COX_SUMMARY_ALL_FILTERED,
        tumor_main_dir = OUTPUT_DIR
    output:
        out_file = FIG_PC_PDL1_EXPRESSION
    message:
        "Generating a figure for the selected cancer types of PDL1 exp vs corresponding PC component"
    params:
        bin = config["bin"],
    shell:
        """
        Rscript {params.bin}/plot_PC_PDL1_exp.R \
            --cox_summary_all_cancers {input.cox_summary_all} \
            --tumor_dir {input.tumor_main_dir} \
            --output {output.out_file}
        """

## Plot the correlation of the PC components and the immune cell types ##

rule plot_PC_immune_correlations:
    input:
        cox_summary_all = UNIVARIATE_COX_SUMMARY_ALL_FILTERED,
        tumor_main_dir = OUTPUT_DIR
    output:
        out_file = FIG_PC_IMMUNE_CORRELATION
    message:
        "Generating a figure of the correlation between the PC and immune cells for selected cancer types"
    params:
        bin = config["bin"],
    shell:
        """
        Rscript {params.bin}/plot_PC_immune_correlations.R \
            --cox_summary_all_cancers {input.cox_summary_all} \
            --tumor_dir {input.tumor_main_dir} \
            --output {output.out_file}
        """
## Perform the association analysis between the PD1 pathway heterogeneity scores and clinical features (categorical and numeric) ##

rule calculate_association_clinical_features_pd1_heterogeneity_scores:
    input:
        cox_res_file = UNIVARIATE_COX_SUMMARY_ALL_FILTERED,
        tumor_main_dir = OUTPUT_DIR
    output:
        results_pd1_groups = TUMOR_RESULTS_PD1_GROUPS,
        results_pd1_numeric = TUMOR_RESULTS_PD1_NUMERIC

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
            --results_pd1_numeric {output.results_pd1_numeric}
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
            --pc_clinical_association_figure {output.figure_pc_clin_associations}

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
            --output_figure_file {output.figure_pc_individual_clin_associations}

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
            --output {output.out_file}
        """
## Combine all multivarite cox results in one table ##  

rule combine_multivarite_cox_results:
    input:
        expand(OUTPUT_CANCER_COX, cancer=CANCER_TYPES)
    output:
        COX_RESULTS_ALL_MULTIVARIATE
    message:
        "Combining all multivariate Cox regression results into one table."
    shell:
        """
        echo -e "cancer\\t$(head -n 1 {input[0]})" > {output}
        for file in {input}; do
            cancer=$(basename $(dirname $(dirname $file)))
            tail -n +2 $file | awk -v cancer=$cancer '{{print cancer"\\t"$0}}' >> {output}
        done
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
            --output {output.out_file}
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
    params:
        bin = config["bin"],
    shell:
        """
        Rscript {params.bin}/plot_TSNE_all_cancers.R \
            --datasets_to_plot_cola_clusters {input.datasets_to_plot_cola_clusters} \
            --tsne_data_expression {input.tsne_file_expression} \
            --tsne_data_indegree {input.tsne_file_indegree} \
            --cancer_color_file {input.cancer_color_file} \
            --ouput_figure_file {output.ouput_figure}
        """
## plot PD1 summary table ##
rule plot_PD1_summary_table:
    input:
        summary_table_PD1 = SUMMARY_TABLE_PD1,
    output:
        output_html_file = OUTPUT_HTML_TABLE_PD1
    params:
        bin = config["bin"],
    shell:
        """
        Rscript {params.bin}/plot_PD1_summary_table.R \
            --summary_table_PD1 {input.summary_table_PD1} \
            --output_html_file {output.output_html_file} 
        """
