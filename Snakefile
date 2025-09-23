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
OUTPUT_DIR = config["results_dir"]
OUTPUT_DIR_INDIVIDUAL_CANCERS = config["output_dir_individual_cancers"]
OUTPUT_DIR_ALL_CANCERS = config["output_dir_all_cancers"]
FIG_DIR = config["fig_dir"]
ENV_DIR = config["envs"]
## From config ##
CANCER_TYPES = config["cancer_types"]
BATCH_FILE = config["batch_file"]
DATATYPES = config["datatypes"]
CLINICAL_FILE = config["clinical_file"]
CLINICAL_FILE_RDATA = config["clinical_file_rdata"]
CANCER_COLOR_FILE = config["cancer_color_file"]
PPI_FILE = config["ppi_file"]
MOTIF_FILE = config["motif_file"]
SAMPLES_FILE = config["samples_file"]
IMMUNE_FILE = config["immune_file"]
GMT_FILE = config["gmt_file"]
PATHWAYS_HIERARCHY_FILE = config["pathways_hierarchy_file"]
PATHWAYS_HSA_ID_FILE = config["pathways_hsa_id_file"]
LIST_PATHWAYS_FILE = config["list_of_pathways_file"]
ANALYSIS_TYPE = config["analysis_type"]
ZENODO_INDIVIDUAL_CANCERS_DIRECTORY = config["zenodo_individual_cancers_directory"]

## Zenodo configuration ##
ZENODO_RECORD_ID = config["zenodo_record_id"]
ZENODO_RESOURCE_FILENAME = config["zenodo_resource_filename"]

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
NCORES_PORCUPINE = config["ncores_porcupine"]


##############################################################################
### ZENODO RESOURCE DOWNLOAD PATHS                                        ###
###############################################################################

## Marker file to indicate successful download of Zenodo resources
ZENODO_RESOURCES_DOWNLOAD_COMPLETE = ".zenodo_download_complete"

ZENODO_BATCH_FILENAME = config["zenodo_batch_analysis_filename"]
## Marker file to indicate successful download of Zenodo resources
ZENODO_BATCH_DOWNLOAD_COMPLETE = ".zenodo_batch_download_complete"


###############################################################################
### DATA DOWNLOAD AND NORMALIZATION PIPELINE PATHS                        ###
###############################################################################

## Primary output directory for downloaded TCGA GDC data
OUTPUT_GDC_DIR = os.path.join(OUTPUT_DIR_ALL_CANCERS, "gdc_data")

## Individual cancer type expression GDC files (RData format)
## Pattern: TCGA-{cancer}.RData (e.g., TCGA-BRCA.RData, TCGA-LUAD.RData)
OUTPUT_GDC_FILE = os.path.join(OUTPUT_GDC_DIR, "TCGA-{cancer}.RData")

#------------------------------------------------------------------------------
# Pre-downloaded Data and Combined Processing
#------------------------------------------------------------------------------
## Directory containing pre-downloaded GDC data (alternative to fresh download)
## Use this when download_files="NO" to work with existing data
EXPRESSION_DIR_GDC = os.path.join(OUTPUT_DIR_ALL_CANCERS, "gdc_data")

## Output directory for combined multi-cancer datasets
OUTPUT_DIR_DOWNLOAD_COMBINED = os.path.join(OUTPUT_DIR_ALL_CANCERS, "combined_gdc_data")

## Combined expression matrix across all cancer types
## Contains raw STAR gene counts in TSV format (genes x samples)
OUTPUT_EXP_COMBINED_FILE = os.path.join(OUTPUT_DIR_DOWNLOAD_COMBINED, "hg38_STAR_counts.tsv")

## Sample-to-cancer type mapping file
## Maps individual sample IDs to their respective cancer types
GROUP_FILE = os.path.join(OUTPUT_DIR_DOWNLOAD_COMBINED, "hg38_sample_groups.tsv")

## Gene annotation and feature information (RData format)
## Contains gene IDs, symbols, biotypes, and genomic coordinates
FEATURE_FILE = os.path.join(OUTPUT_DIR_DOWNLOAD_COMBINED, "hg38_features.RData")

#------------------------------------------------------------------------------
# Normalized Expression Data
#------------------------------------------------------------------------------
## PySNAIL normalized expression file (all cancers combined)
## Contains qsmooth-normalized gene expression values
PYSNAIL_NORMALIZED_FILE = os.path.join(OUTPUT_DIR_DOWNLOAD_COMBINED, "pysnail_normalized_STAR_counts.tsv")

## Directory for cancer-specific normalized expression files
## Contains individual RData files for each cancer type post-normalization
## Pattern: normalized_expression_TCGA-{cancer}.RData
OUTPUT_DIR_PYSNAIL_CANCER = os.path.join(OUTPUT_DIR_ALL_CANCERS, "pysnail_normalized_individual_cancer_expression")
PYSNAIL_NORMALIZED_FILE_CANCER_SPECIFIC = os.path.join(OUTPUT_DIR_PYSNAIL_CANCER, "normalized_expression_TCGA-{cancer}.RData")
PYSNAIL_YAML_RELATIVE = os.path.relpath(
    os.path.join("workflow", "envs", "pysnail.yaml"), 
    "workflow/rules")

###############################################################################
### BATCH EFFECT ANALYSIS PIPELINE PATHS                                  ###
###############################################################################

#------------------------------------------------------------------------------
# Batch Effect Detection and Analysis
#------------------------------------------------------------------------------
## Primary output directory for batch effect analysis results across all cancers
## Contains MBatch DSC analysis results and diagnostic files
BATCH_DIR_ALL_CANCERS = os.path.join(OUTPUT_DIR_ALL_CANCERS, "batch_analysis")

## Cancer-specific batch analysis directories  
## Each directory contains PCA plots, DSC statistics, and batch effect diagnostics
BATCH_DIR_CANCER = os.path.join(BATCH_DIR_ALL_CANCERS, "TCGA-{cancer}")

## Summary figure showing Dispersive Separation Criterion (DSC) values
BATCH_EFFECT_PDF = os.path.join(FIG_DIR, "MBatch_DSC.pdf")

#------------------------------------------------------------------------------
# Batch-Corrected Expression Data
#------------------------------------------------------------------------------
## ComBat-corrected expression matrix (all cancers combined)
## Contains batch-effect corrected log2-transformed gene expression values
BATCH_CORRECTED_EXPRESSION_FILE = os.path.join(
    OUTPUT_DIR_ALL_CANCERS, 
    "batch_corrected_expression", 
    "batch_corrected_expression_all_cancers.RData"
)


#------------------------------------------------------------------------------
# PANDA Prior Networks and Expression Data
#------------------------------------------------------------------------------
## Directory: panda_input/ contains all filtered input files for PANDA/LIONESS

NETZOOPY_YAML = os.path.relpath(
    os.path.join("workflow", "envs", "netzoopy-local.yaml"), 
    "workflow/rules")

## Transcription factor binding motif prior network (filtered)
MOTIF_PANDA_FILE = os.path.join(OUTPUT_DIR, "panda_input/motif_tcga_primary.tsv")

## Protein-protein interaction prior network (filtered)  
PPI_PANDA_FILE = os.path.join(OUTPUT_DIR, "panda_input/ppi_tcga_primary.tsv")

## Gene expression matrix (batch-corrected and filtered)
EXPRESSION_PANDA_FILE = os.path.join(OUTPUT_DIR, "panda_input/exp_tcga_primary.tsv")

## Sample identifier list (primary tumors only)
## Contains TCGA sample IDs in order matching expression matrix columns
SAMPLES_PANDA_FILE = os.path.join(OUTPUT_DIR, "panda_input/samples_primary.tsv")

## Sample-to-cancer mapping file (primary tumors only)
## Maps each sample to its corresponding cancer type for downstream analysis
SAMPLES_WITH_CANCER_FILE = os.path.join(OUTPUT_DIR, "panda_input/samples_cancers_primary.tsv")

###############################################################################
### PANDA + LIONESS NETWORK INFERENCE AND POST-PROCESSING                    ###
###############################################################################

#------------------------------------------------------------------------------
# Raw LIONESS Network Inference Output
#------------------------------------------------------------------------------
## Primary output directory for PANDA aggregate and LIONESS sample-specific networks
## Contains individual .txt files for each sample's regulatory network
NETWORKS_DIR = os.path.join(OUTPUT_DIR_ALL_CANCERS, "networks")

## Comprehensive mapping file linking LIONESS network files to sample metadata
LIONESS_SAMPLE_MAPPING = os.path.join(NETWORKS_DIR, "information_networks_primary.txt")

#------------------------------------------------------------------------------
# Cancer-Specific Network Aggregation
#------------------------------------------------------------------------------
## Directory for cancer-specific merged network RData files
## Contains combined networks for each cancer type in efficient RData format
## Files: net_BRCA.RData, net_LUAD.RData, etc.
## Large cancer types (≥300 samples) are split into multiple files (e.g., net_BRCA1.RData, net_BRCA2.RData)
OUTPUT_DIR_FINAL_MERGED_NETWORKS = os.path.join(OUTPUT_DIR_ALL_CANCERS, "final_networks")

#------------------------------------------------------------------------------
# Quantile-Normalized Networks
#------------------------------------------------------------------------------
## Directory for quantile-normalized cancer-specific networks
## Contains normalized versions of merged networks 
## Files: net_norm_TCGA-BRCA.RData, net_norm_TCGA-LUAD.RData, etc.
OUTPUT_DIR_NORMALIZED_NETWORKS = os.path.join(OUTPUT_DIR_ALL_CANCERS, "final_networks_normalized")

###############################################################################
PANDA_NETWORK_FILE = os.path.join(NETWORKS_DIR, "panda_net.txt") 
NETWORK_EDGE_FILE = os.path.join(OUTPUT_DIR_ALL_CANCERS, "edges", "network_edges.txt")

#------------------------------------------------------------------------------
# Gene Indegree Calculation from Normalized Networks
#------------------------------------------------------------------------------
CANCER_INDEGREE_FILE  = os.path.join(OUTPUT_DIR_INDIVIDUAL_CANCERS, "{cancer}", "indegrees_norm", "indegree_norm_{cancer}.RData")

###############################################################################
### TSNE AND PORCUPINE ANALYSIS OUTPUTS                                     ###
###############################################################################
OUTPUT_CANCER = os.path.join(OUTPUT_DIR_INDIVIDUAL_CANCERS, "{cancer}", "clinical", "curated_clinical_{cancer}.txt")

# TSNE Outputs
# Directory for TSNE results across all cancers
TSNE_DATA_DIR = os.path.join(OUTPUT_DIR_ALL_CANCERS, "tsne_results")

# TSNE results for all cancers based on expression data
TSNE_DATA_EXPRESSION = os.path.join(TSNE_DATA_DIR, "tsne_expression_all_cancers.txt")

# TSNE results for all cancers based on indegree data
TSNE_DATA_INDEGREE = os.path.join(TSNE_DATA_DIR, "tsne_expression_all_indegree.txt")


# Filtered Porcupine results for each cancer type
# Porcupine results with variance for each cancer type
PORCUPINE_FILE = os.path.join(OUTPUT_DIR_INDIVIDUAL_CANCERS, "{cancer}", "porcupine", "pcp_results_with_variance_{cancer}.txt")

FILTERED_PORCUPINE_FILE = os.path.join(OUTPUT_DIR_INDIVIDUAL_CANCERS, "{cancer}", "porcupine", "pcp_results_filtered_{cancer}.txt")

# Combined Porcupine results for all cancers
PORCUPINE_RESULTS_ALL = os.path.join(OUTPUT_DIR_ALL_CANCERS, "porcupine", "porcupine_results_all_cancers.txt")

# PDF figure showing the intersection of pathways across cancers
FIG_PATHWAY_INTERSECTION = os.path.join(FIG_DIR, "pathways_intersection_pcp.pdf")

# PDF figure showing shared categories across cancers
FIG_SHARED_CATEGORIES = os.path.join(FIG_DIR, "combined_figure_pcp_results_shared_categories.pdf")

# Curated clinical data for each cancer type
TUMOR_CLIN_FILE = os.path.join(OUTPUT_DIR_INDIVIDUAL_CANCERS, "{cancer}", "clinical", "curated_clinical_{cancer}.txt")

# PDF figure showing the legend for cancer types
CANCER_LEGEND_PDF = os.path.join(FIG_DIR, "cancer_legend.pdf")

# PDL1 expression data for each cancer type
INPUT_CANCER_INDEGREE_DIR  = os.path.join(OUTPUT_DIR_INDIVIDUAL_CANCERS, "{cancer}", "indegrees_norm")
TUMOR_PD1_DIR = os.path.join(OUTPUT_DIR_INDIVIDUAL_CANCERS, "{cancer}", "pd1_data")
OUTPUT_PDL1_EXP_CANCER = os.path.join(OUTPUT_DIR_INDIVIDUAL_CANCERS, "{cancer}", "pd1_data", "pdl1_expression_{cancer}.txt")
OUTPUT_CANCER_PD1_MAPPINGS  = os.path.join(OUTPUT_DIR_INDIVIDUAL_CANCERS, "{cancer}", "pd1_data", "pd1_individual_scores_norm_{cancer}.RData")
TUMOR_PATHWAYS_MAPPING_PATH = os.path.join(OUTPUT_DIR_INDIVIDUAL_CANCERS, "{cancer}", "porcupine", "individual_scores_{cancer}.RData")


###############################################################################
### COLA CONSENSUS CLUSTERING ANALYSIS PIPELINE PATHS                     ###
###############################################################################
#------------------------------------------------------------------------------
# Cancer-Specific Consensus Clustering Results
#------------------------------------------------------------------------------
## Individual cancer type consensus clustering output directories
## Pattern: {cancer}/consensus_clustering/{datatype} (e.g., BRCA/consensus_clustering/expression)
## Contains COLA analysis results for each cancer type and data type combination
OUTPUT_CANCER_CONSENSUS_DIR = os.path.join(OUTPUT_DIR_INDIVIDUAL_CANCERS, "{cancer}", "consensus_clustering", "{datatype}")
BEST_K_COLA_OUTPUT_CANCER_DATATYPE = os.path.join(OUTPUT_CANCER_CONSENSUS_DIR, "{cancer}_best_k_{datatype}.RData")
COLA_RESULTS_OUTPUT_CANCER_DATATYPE = os.path.join(OUTPUT_CANCER_CONSENSUS_DIR, "{cancer}_results_{datatype}.RData")
COLA_COLLECTED_PLOTS_CANCER_DATATYPE = os.path.join(OUTPUT_CANCER_CONSENSUS_DIR, "{cancer}_collected_plots_{datatype}.pdf")
COLA_TSNE_OUTPUT_CANCER_DATATYPE = os.path.join(OUTPUT_CANCER_CONSENSUS_DIR, "{cancer}_tSNE_{datatype}.pdf")     
COLA_MEMBERSHIP_OUTPUT_CANCER_DATATYPE = os.path.join(OUTPUT_CANCER_CONSENSUS_DIR, "{cancer}_membership_{datatype}.RData")
COLA_STATISTICS_OUTPUT_CANCER_DATATYPE = os.path.join(OUTPUT_CANCER_CONSENSUS_DIR, "{cancer}_statistics_{datatype}.RData")
COLA_PARTITION_OUTPUT_CANCER_DATATYPE = os.path.join(OUTPUT_CANCER_CONSENSUS_DIR, "{cancer}_select_partition_{datatype}.pdf")
COLA_CLASSES_OUTPUT_CANCER_DATATYPE = os.path.join(OUTPUT_CANCER_CONSENSUS_DIR, "{cancer}_classes_{datatype}.RData")

OUTPUT_CANCER_CONSENSUS_DIR_INDEGREE = os.path.join(OUTPUT_DIR_INDIVIDUAL_CANCERS, "{cancer}", "consensus_clustering", "indegree")
OUTPUT_CANCER_CONSENSUS_DIR_EXPRESSION = os.path.join(OUTPUT_DIR_INDIVIDUAL_CANCERS, "{cancer}", "consensus_clustering", "expression")

BEST_K_COLA_INDGEGREE_FILES = expand(
    os.path.join(OUTPUT_CANCER_CONSENSUS_DIR_INDEGREE, "{cancer}_best_k_indegree.RData"),
    cancer=CANCER_TYPES
)

BEST_K_COLA_EXPRESSION_FILES = expand(
    os.path.join(OUTPUT_CANCER_CONSENSUS_DIR_EXPRESSION, "{cancer}_best_k_expression.RData"),
    cancer=CANCER_TYPES
)

COLA_RESULTS_INDEGREE_FILES = expand(
    os.path.join(OUTPUT_CANCER_CONSENSUS_DIR_INDEGREE, "{cancer}_results_indegree.RData"),
    cancer=CANCER_TYPES
)

COLA_RESULTS_EXPRESSION_FILES = expand(
    os.path.join(OUTPUT_CANCER_CONSENSUS_DIR_EXPRESSION, "{cancer}_results_expression.RData"),
    cancer=CANCER_TYPES
)

#------------------------------------------------------------------------------
# Pan-Cancer Consensus Clustering Analysis
#------------------------------------------------------------------------------
## Primary output directory for combined COLA consensus clustering results
## Contains optimized cluster assignments and comparative analysis across all cancer types
OUTPUT_ALL_CANCERS_CONSENSUS_DIR = os.path.join(OUTPUT_DIR_ALL_CANCERS, "cola_consensus_clustering")

## Optimal cluster number determination files
## Contains best K values selected by COLA algorithm for each data type
BEST_K_COLA_EXP = os.path.join(OUTPUT_ALL_CANCERS_CONSENSUS_DIR, "best_k_cola_expression.txt")
BEST_K_COLA_IND = os.path.join(OUTPUT_ALL_CANCERS_CONSENSUS_DIR, "best_k_cola_indegree.txt")

## Sample-to-cluster assignment tables for selected optimal K values
## Maps individual samples to their assigned cluster groups for downstream analysis
SELECTED_CLUSTERS_COLA_EXP = os.path.join(OUTPUT_ALL_CANCERS_CONSENSUS_DIR, "selected_clusters_expression.txt")
SELECTED_CLUSTERS_COLA_IND = os.path.join(OUTPUT_ALL_CANCERS_CONSENSUS_DIR, "selected_clusters_indegree.txt")

## Combined dataset for cluster visualization and comparative analysis
## RData object containing processed data for t-SNE and other visualization methods
DATASETS_TO_PLOT_COLA_CLUSTERS = os.path.join(OUTPUT_ALL_CANCERS_CONSENSUS_DIR, "datasets_to_plot_cola_clusters.RData")

#------------------------------------------------------------------------------
# Consensus Clustering Visualization Outputs
#------------------------------------------------------------------------------
## t-SNE visualization of consensus clusters overlaid on dimensional reduction plots
## Shows sample clustering patterns and cluster separation quality
FIG_TSNE_COLA_INDEGREE = os.path.join(FIG_DIR, "TSNE_cola_clusters_indegree_all_cancers.pdf")
FIG_TSNE_COLA_EXPRESSION = os.path.join(FIG_DIR, "TSNE_cola_clusters_expression_all_cancers.pdf")

## Sankey diagram comparing cluster assignments between expression and indegree data
## Visualizes concordance and differences between clustering approaches
FIG_SANKEY = os.path.join(FIG_DIR, "sankey_plot_indegree_expression.pdf")

#------------------------------------------------------------------------------
# Final Cancer-Specific Cluster Assignments
#------------------------------------------------------------------------------
## Final cluster assignment files for individual cancer types
## Contains optimized cluster memberships for downstream survival and clinical analysis
OUTPUT_CLUSTERS_PER_TUMOR_IND = os.path.join(OUTPUT_DIR_INDIVIDUAL_CANCERS, "{cancer}", "final_clusters", "final_clusters_indegree_{cancer}.txt")
OUTPUT_CLUSTERS_PER_TUMOR_EXP = os.path.join(OUTPUT_DIR_INDIVIDUAL_CANCERS, "{cancer}", "final_clusters", "final_clusters_expression_{cancer}.txt")

#------------------------------------------------------------------------------
# Survival Analysis of Consensus Clusters
#------------------------------------------------------------------------------
## Cox regression results comparing survival outcomes between consensus clusters
## Individual cancer type analysis results
OUTPUT_CANCER_UNIVARIATE_COX_COLA_CLUSTERS = os.path.join(OUTPUT_DIR_INDIVIDUAL_CANCERS, "{cancer}", "final_clusters", "cox_results_final_clusters_indegree_expression_{cancer}.txt")

## Combined Cox regression results across all cancer types
## Comprehensive survival analysis summary for consensus cluster validation
OUTPUT_CANCER_UNIVARIATE_COX_COLA_CLUSTERS_ALL = os.path.join(OUTPUT_DIR_ALL_CANCERS, "cox_results_all", "cox_results_final_clusters_indegree_expression_all.txt")

## Visualization of Cox regression results across cancer types
## Summary figure showing prognostic significance of consensus clusters
FIG_COX_COLA_CLUSTERS = os.path.join(FIG_DIR, "cox_results_final_clusters_indegree_expression_all_cancers.pdf") 

#------------------------------------------------------------------------------
# Aanalysis of PRAD Clusters at K = 4
#------------------------------------------------------------------------------

CLUSTER_INDEGREE_PRAD = os.path.join(OUTPUT_DIR_INDIVIDUAL_CANCERS, "PRAD", "final_clusters", "final_clusters_indegree_PRAD.txt")
PRAD_CLIN_FILE = os.path.join(OUTPUT_DIR_INDIVIDUAL_CANCERS, "PRAD", "clinical", "curated_clinical_PRAD.txt")
PRAD_PD1_DIR = os.path.join(OUTPUT_DIR_INDIVIDUAL_CANCERS, "PRAD", "pd1_data")
PRAD_IND_FILE = os.path.join(OUTPUT_DIR_INDIVIDUAL_CANCERS, "PRAD", "indegrees_norm", "indegree_norm_PRAD.RData")
FIG_PRAD_SURVIVAL = os.path.join(FIG_DIR, "PRAD_clusters_survival.pdf")
FIG_FGSEA_PRAD = os.path.join(FIG_DIR, "PRAD_clusters_fgsea.pdf")


## Output Files for Univariate Cox on PD1-pathway based heterogeneity scores ##
OUTPUT_CANCER_UNIVARIATE_COX_SUMMARY = os.path.join(OUTPUT_DIR_INDIVIDUAL_CANCERS, "{cancer}", "cox", "{cancer}_PD1_pathway_cox_univariate_model_summary.txt")
OUTPUT_CANCER_UNIVARIATE_COX_PREDICTED_SCORES = os.path.join(OUTPUT_DIR_INDIVIDUAL_CANCERS, "{cancer}", "cox", "{cancer}_PD1_pathway_cox_univariate_predited_risk_scores.txt")
UNIVARIATE_COX_SUMMARY_ALL = os.path.join(OUTPUT_DIR_ALL_CANCERS, "cox_results_all", "PD1_pathway_cox_univariate_model_summary_all.txt")
UNIVARIATE_COX_SUMMARY_ALL_FILTERED = os.path.join(OUTPUT_DIR_ALL_CANCERS, "cox_results_all", "PD1_pathway_cox_univariate_model_summary_filtered.txt")
UNIVARIATE_COX_PREDICTED_SCORES_ALL = os.path.join(OUTPUT_DIR_ALL_CANCERS, "cox_results_all", "PD1_pathway_cox_univariate_predited_risk_scores_all.txt")
FIG_PC_PDL1_EXPRESSION = os.path.join(FIG_DIR, "PDL1_exp_PC_component_HR.pdf")
FIG_PC_IMMUNE_CORRELATION = os.path.join(FIG_DIR, "PC_immune_correlations_cibersort.png")
TUMOR_RESULTS_PD1_GROUPS = os.path.join(OUTPUT_DIR_ALL_CANCERS, "clinical_associations_PD1", "pd1_pathway_categorical_results.txt")
TUMOR_RESULTS_PD1_NUMERIC = os.path.join(OUTPUT_DIR_ALL_CANCERS, "clinical_associations_PD1", "pd1_pathway_numeric_results.txt")
FIGURE_PC_CLIN_ASSOCIATIONS = os.path.join(FIG_DIR, "PC_all_features_clin_associations.pdf")
FIGURE_PC_INDIVIDUAL_CLIN_ASSOCIATIONS = os.path.join(FIG_DIR, "PC_individual_features_clin_associations.pdf")

## Output Files for multivariate regularized Cox on PDL1-edges ##
OUTPUT_CANCER_PD1_MAPPINGS  = os.path.join(OUTPUT_DIR_INDIVIDUAL_CANCERS, "{cancer}", "pd1_data", "pd1_individual_scores_norm_{cancer}.RData")
OUTPUT_CANCER_COX = os.path.join(OUTPUT_DIR_INDIVIDUAL_CANCERS, "{cancer}", "cox", "{cancer}_PDL1_cox_multivariate_res.txt")
COX_RESULTS_ALL_MULTIVARIATE = os.path.join(OUTPUT_DIR_ALL_CANCERS, "cox_results_all", "PDL1_cox_multivarite_res_all.txt")
PDL1_CIRCULAR_PLOT = os.path.join(FIG_DIR, "circular_pdl1_plot_{threshold_cox}.pdf")
OUTPUT_PDL1_EXP_CANCER = os.path.join(OUTPUT_DIR_INDIVIDUAL_CANCERS, "{cancer}", "pd1_data", "pdl1_expression_{cancer}.txt")
OUTPUT_COMBINED_PATIENT_DATA_CANCER = os.path.join(OUTPUT_DIR_INDIVIDUAL_CANCERS, "{cancer}", "pd1_data", "combined_patient_data_{cancer}.txt")


# figure for TNSE plot for all cancers and also comparisons of cola clusters for indegree and expression for PRAD and UVM
FIGURE_TSNE_ALL_CANCERS_UVM_PRAD_CLUSTERS = os.path.join(FIG_DIR, "TSNE_all_cancers_indegree_expression_UVM_PRAD_clusters.pdf")

# producing summary table figure for PD1 pathway 
#input (manually created)
SUMMARY_TABLE_PD1 = os.path.join("resources", "summary_table", "summary_table_PD1.txt")
#output
OUTPUT_HTML_TABLE_PD1 = os.path.join(FIG_DIR, "summary_table_PD1.html")




# Always include all rule files (Snakemake will only execute needed rules)
include: "workflow/rules/zenodo_download_resources.smk"
include: "workflow/rules/batch_effect_expression.smk"
include: "workflow/rules/download_normalize_expression_data.smk"
include: "workflow/rules/zenodo_download_batch_results.smk"
include: "workflow/rules/prepare_data_run_networks.smk"
include: "workflow/rules/zenodo_download_indegrees_porcupine_data.smk"
include: "workflow/rules/tsne_and_porcupine_analysis.smk"
include: "workflow/rules/extract_pd1_data.smk"
include: "workflow/rules/cola_consensus_clustering.smk"
include: "workflow/rules/prad_cluster_analysis.smk"
include: "workflow/rules/pd1_analysis.smk"

# Rules ##
rule all:
    input:
        # Always download resources (both workflows need this)
        # ZENODO_RESOURCES_DOWNLOAD_COMPLETE,      
        # Conditional inputs based on analysis type
        # ([ZENODO_BATCH_DOWNLOAD_COMPLETE] if ANALYSIS_TYPE == "precomputed" else []),
        # Full workflow specific outputs
        ([expand(OUTPUT_GDC_FILE, cancer=CANCER_TYPES)] if ANALYSIS_TYPE == "full_workflow" else []),
        ([OUTPUT_EXP_COMBINED_FILE] if ANALYSIS_TYPE == "full_workflow" else []),
        ([GROUP_FILE] if ANALYSIS_TYPE == "full_workflow" else []),
        ([FEATURE_FILE] if ANALYSIS_TYPE == "full_workflow" else []),
        ([PYSNAIL_NORMALIZED_FILE] if ANALYSIS_TYPE == "full_workflow" else []),
        ([expand(PYSNAIL_NORMALIZED_FILE_CANCER_SPECIFIC, cancer = CANCER_TYPES)] if ANALYSIS_TYPE == "full_workflow" else []),
        ([expand(BATCH_DIR_CANCER, cancer=CANCER_TYPES)] if ANALYSIS_TYPE == "full_workflow" else []),
        # ([BATCH_CORRECTED_EXPRESSION_FILE] if ANALYSIS_TYPE == "full_workflow" else []),       
        # # # Shared outputs (both workflows can generate this)
        # ([BATCH_EFFECT_PDF] if ANALYSIS_TYPE in ["full_workflow", "precomputed"] else []),
        # # PANDA/LIONESS files (only for full workflow),
        ([MOTIF_PANDA_FILE] if ANALYSIS_TYPE in ["full_workflow"] else []),
        ([PPI_PANDA_FILE] if ANALYSIS_TYPE in ["full_workflow"] else []),
        ([EXPRESSION_PANDA_FILE] if ANALYSIS_TYPE in ["full_workflow"] else []),
        ([SAMPLES_PANDA_FILE] if ANALYSIS_TYPE in ["full_workflow"] else []),
        ([SAMPLES_WITH_CANCER_FILE] if ANALYSIS_TYPE in ["full_workflow"] else []),
        ## Download indegree files and porcupine data from Zenodo 
        # ([expand(ZENODO_INDIVIDUAL_CANCERS_DIRECTORY)] if ANALYSIS_TYPE in ["precomputed"] else []) 
        # HERE IS PART TO RUN ONLY FOR FULL WORKFLOW
        # ([NETWORKS_DIR] if ANALYSIS_TYPE in ["full_workflow"] else []),
        # ([LIONESS_SAMPLE_MAPPING] if ANALYSIS_TYPE in ["full_workflow"] else []),
        # ([OUTPUT_DIR_FINAL_MERGED_NETWORKS] if ANALYSIS_TYPE in ["full_workflow"] else []),
        # ([OUTPUT_DIR_NORMALIZED_NETWORKS] if ANALYSIS_TYPE in ["full_workflow"] else []),
        # ([PANDA_NETWORK_FILE] if ANALYSIS_TYPE in ["full_workflow"] else []),
        # ([NETWORK_EDGE_FILE] if ANALYSIS_TYPE in ["full_workflow"] else []),
        # ([expand(CANCER_INDEGREE_FILE, cancer=CANCER_TYPES)] if ANALYSIS_TYPE in ["full_workflow"] else [])
         # HERE IS PART TO RUN INDEPENDENTLY OF ANALYSIS TYPE
        ([expand(OUTPUT_CANCER, cancer=CANCER_TYPES)] if ANALYSIS_TYPE in ["full_workflow", "precomputed"] else []),
        # TSNE Outputs
        ([TSNE_DATA_EXPRESSION] if ANALYSIS_TYPE in ["full_workflow", "precomputed"] else []),
        ([TSNE_DATA_INDEGREE] if ANALYSIS_TYPE in ["full_workflow", "precomputed"] else []),
        ([CANCER_LEGEND_PDF] if ANALYSIS_TYPE in ["full_workflow", "precomputed"] else []),
        # Porcupine Outputs
        ([expand(FILTERED_PORCUPINE_FILE, cancer = CANCER_TYPES)] if ANALYSIS_TYPE in ["full_workflow", "precomputed"] else []),
        ([PORCUPINE_RESULTS_ALL] if ANALYSIS_TYPE in ["full_workflow", "precomputed"] else []),
        ([FIG_PATHWAY_INTERSECTION] if ANALYSIS_TYPE in ["full_workflow", "precomputed"] else []),
        ([FIG_SHARED_CATEGORIES] if ANALYSIS_TYPE in ["full_workflow", "precomputed"] else []),
        ([expand(TUMOR_CLIN_FILE, cancer = CANCER_TYPES)] if ANALYSIS_TYPE in ["full_workflow", "precomputed"] else []),
        # cola consensus clustering outputs
        ([expand(BEST_K_COLA_OUTPUT_CANCER_DATATYPE, cancer = CANCER_TYPES, datatype = DATATYPES)] if ANALYSIS_TYPE in ["full_workflow", "precomputed"] else []),
        ([expand(COLA_RESULTS_OUTPUT_CANCER_DATATYPE, cancer = CANCER_TYPES, datatype = DATATYPES)] if ANALYSIS_TYPE in ["full_workflow", "precomputed"] else []),
        ([expand(COLA_COLLECTED_PLOTS_CANCER_DATATYPE, cancer = CANCER_TYPES, datatype = DATATYPES)] if ANALYSIS_TYPE in ["full_workflow", "precomputed"] else []),
        ([expand(COLA_TSNE_OUTPUT_CANCER_DATATYPE, cancer = CANCER_TYPES, datatype = DATATYPES)] if ANALYSIS_TYPE in ["full_workflow", "precomputed"] else []),
        ([expand(COLA_MEMBERSHIP_OUTPUT_CANCER_DATATYPE, cancer = CANCER_TYPES, datatype = DATATYPES)] if ANALYSIS_TYPE in ["full_workflow", "precomputed"] else []),
        ([expand(COLA_STATISTICS_OUTPUT_CANCER_DATATYPE, cancer = CANCER_TYPES, datatype = DATATYPES)] if ANALYSIS_TYPE in ["full_workflow", "precomputed"] else []),
        ([expand(COLA_PARTITION_OUTPUT_CANCER_DATATYPE, cancer = CANCER_TYPES, datatype = DATATYPES)] if ANALYSIS_TYPE in ["full_workflow", "precomputed"] else []),
        ([expand(COLA_CLASSES_OUTPUT_CANCER_DATATYPE, cancer = CANCER_TYPES, datatype = DATATYPES)] if ANALYSIS_TYPE in ["full_workflow", "precomputed"] else []),
        ([BEST_K_COLA_EXP] if ANALYSIS_TYPE in ["full_workflow", "precomputed"] else []),
        ([BEST_K_COLA_IND] if ANALYSIS_TYPE in ["full_workflow", "precomputed"] else []),
        ([FIG_TSNE_COLA_INDEGREE] if ANALYSIS_TYPE in ["full_workflow", "precomputed"] else []),
        ([FIG_TSNE_COLA_EXPRESSION] if ANALYSIS_TYPE in ["full_workflow", "precomputed"] else []),
        ([FIG_SANKEY] if ANALYSIS_TYPE in ["full_workflow", "precomputed"] else []),
        ([SELECTED_CLUSTERS_COLA_EXP] if ANALYSIS_TYPE in ["full_workflow", "precomputed"] else []),
        ([SELECTED_CLUSTERS_COLA_IND] if ANALYSIS_TYPE in ["full_workflow", "precomputed"] else []),
        ([expand(OUTPUT_CLUSTERS_PER_TUMOR_IND, cancer = CANCER_TYPES)] if ANALYSIS_TYPE in ["full_workflow", "precomputed"] else []),
        ([expand(OUTPUT_CLUSTERS_PER_TUMOR_EXP, cancer = CANCER_TYPES)] if ANALYSIS_TYPE in ["full_workflow", "precomputed"] else []),
        # cox analysis of cola clusters
        ([expand(OUTPUT_PDL1_EXP_CANCER, cancer = CANCER_TYPES)] if ANALYSIS_TYPE in ["full_workflow", "precomputed"] else []),
        ([expand(OUTPUT_CANCER_PD1_MAPPINGS, cancer = CANCER_TYPES)] if ANALYSIS_TYPE in ["full_workflow", "precomputed"] else []),
        ([expand(OUTPUT_CANCER_UNIVARIATE_COX_COLA_CLUSTERS, cancer = CANCER_TYPES)] if ANALYSIS_TYPE in ["full_workflow", "precomputed"] else []),
        ([OUTPUT_CANCER_UNIVARIATE_COX_COLA_CLUSTERS_ALL] if ANALYSIS_TYPE in ["full_workflow", "precomputed"] else []),
        ([FIG_COX_COLA_CLUSTERS] if ANALYSIS_TYPE in ["full_workflow", "precomputed"] else []),
        # # prad cluster analysis
        ([FIG_PRAD_SURVIVAL] if ANALYSIS_TYPE in ["full_workflow", "precomputed"] else []),
        ([FIG_FGSEA_PRAD] if ANALYSIS_TYPE in ["full_workflow", "precomputed"] else []),
        # pd1 analysis             
        ([expand(OUTPUT_CANCER_UNIVARIATE_COX_SUMMARY, cancer = CANCER_TYPES)] if ANALYSIS_TYPE in ["full_workflow", "precomputed"] else []),
        ([expand(OUTPUT_CANCER_UNIVARIATE_COX_PREDICTED_SCORES, cancer = CANCER_TYPES)] if ANALYSIS_TYPE in ["full_workflow", "precomputed"] else []),
        ([UNIVARIATE_COX_SUMMARY_ALL] if ANALYSIS_TYPE in ["full_workflow", "precomputed"] else []),    
        ([UNIVARIATE_COX_PREDICTED_SCORES_ALL] if ANALYSIS_TYPE in ["full_workflow", "precomputed"] else []),
        ([UNIVARIATE_COX_SUMMARY_ALL_FILTERED] if ANALYSIS_TYPE in ["full_workflow", "precomputed"] else []),
        ([expand(OUTPUT_COMBINED_PATIENT_DATA_CANCER, cancer = CANCER_TYPES)] if ANALYSIS_TYPE in ["full_workflow", "precomputed"] else []),
        ([FIG_PC_PDL1_EXPRESSION] if ANALYSIS_TYPE in ["full_workflow", "precomputed"] else []),
        ([FIG_PC_IMMUNE_CORRELATION] if ANALYSIS_TYPE in ["full_workflow", "precomputed"] else []),
        ([TUMOR_RESULTS_PD1_GROUPS] if ANALYSIS_TYPE in ["full_workflow", "precomputed"] else []),
        ([TUMOR_RESULTS_PD1_NUMERIC] if ANALYSIS_TYPE in ["full_workflow", "precomputed"] else []),
        ([FIGURE_PC_CLIN_ASSOCIATIONS] if ANALYSIS_TYPE in ["full_workflow", "precomputed"] else []),
        ([FIGURE_PC_INDIVIDUAL_CLIN_ASSOCIATIONS] if ANALYSIS_TYPE in ["full_workflow", "precomputed"] else []),
        ([expand(OUTPUT_CANCER_COX, cancer = CANCER_TYPES)] if ANALYSIS_TYPE in ["full_workflow", "precomputed"] else []),
        ([COX_RESULTS_ALL_MULTIVARIATE] if ANALYSIS_TYPE in ["full_workflow", "precomputed"] else []), 
        ([expand(PDL1_CIRCULAR_PLOT, threshold_cox = THESHOLD_COX)] if ANALYSIS_TYPE in ["full_workflow", "precomputed"] else []),
        ([FIGURE_TSNE_ALL_CANCERS_UVM_PRAD_CLUSTERS] if ANALYSIS_TYPE in ["full_workflow", "precomputed"] else []),
        ([OUTPUT_HTML_TABLE_PD1] if ANALYSIS_TYPE in ["full_workflow", "precomputed"] else [])
        