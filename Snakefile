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


# TUMOR_PD1_DIR = os.path.join(OUTPUT_DIR_INDIVIDUAL_CANCERS, "{cancer}", "pd1_data")
# TUMOR_PATHWAYS_MAPPING_PATH = os.path.join(OUTPUT_DIR_INDIVIDUAL_CANCERS, "{cancer}", "porcupine", "individual_scores_{cancer}.RData")
# INPUT_CANCER_INDEGREE_DIR  = os.path.join(OUTPUT_DIR_INDIVIDUAL_CANCERS, "{cancer}", "indegrees_norm")
# TSNE_DIR = os.path.join(OUTPUT_DIR_ALL_CANCERS, "tsne_results")



# Always include all rule files (Snakemake will only execute needed rules)
include: "workflow/rules/zenodo_download_resources.smk"
include: "workflow/rules/batch_effect_expression.smk"
include: "workflow/rules/download_normalize_expression_data.smk"
include: "workflow/rules/zenodo_download_batch_results.smk"
include: "workflow/rules/prepare_data_run_networks.smk"
include: "workflow/rules/zenodo_download_indegrees_porcupine_data.smk"
include: "workflow/rules/tsne_and_porcupine_analysis.smk"

# Rules ##
rule all:
    input:
        # Always download resources (both workflows need this)
        # ZENODO_RESOURCES_DOWNLOAD_COMPLETE,      
        # Conditional inputs based on analysis type
        # ([ZENODO_BATCH_DOWNLOAD_COMPLETE] if ANALYSIS_TYPE == "precomputed" else []),
        # Full workflow specific outputs
        # ([expand(OUTPUT_GDC_FILE, cancer=CANCER_TYPES)] if ANALYSIS_TYPE == "full_workflow" else []),
        # ([OUTPUT_EXP_COMBINED_FILE] if ANALYSIS_TYPE == "full_workflow" else []),
        # ([GROUP_FILE] if ANALYSIS_TYPE == "full_workflow" else []),
        # ([FEATURE_FILE] if ANALYSIS_TYPE == "full_workflow" else []),
        # ([PYSNAIL_NORMALIZED_FILE] if ANALYSIS_TYPE == "full_workflow" else []),
        # ([OUTPUT_DIR_PYSNAIL_CANCER] if ANALYSIS_TYPE == "full_workflow" else []),
        # ([expand(BATCH_DIR_CANCER, cancer=CANCER_TYPES)] if ANALYSIS_TYPE == "full_workflow" else []),
        # ([BATCH_CORRECTED_EXPRESSION_FILE] if ANALYSIS_TYPE == "full_workflow" else []),       
        # # # Shared outputs (both workflows can generate this)
        # ([BATCH_EFFECT_PDF] if ANALYSIS_TYPE in ["full_workflow", "precomputed"] else []),
        # # PANDA/LIONESS files (only for full workflow),
        # ([MOTIF_PANDA_FILE] if ANALYSIS_TYPE in ["full_workflow"] else []),
        # ([PPI_PANDA_FILE] if ANALYSIS_TYPE in ["full_workflow"] else []),
        # ([EXPRESSION_PANDA_FILE] if ANALYSIS_TYPE in ["full_workflow"] else []),
        # ([SAMPLES_PANDA_FILE] if ANALYSIS_TYPE in ["full_workflow"] else []),
        # ([SAMPLES_WITH_CANCER_FILE] if ANALYSIS_TYPE in ["full_workflow"] else []),
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
        ([expand(TUMOR_CLIN_FILE, cancer = CANCER_TYPES)] if ANALYSIS_TYPE in ["full_workflow", "precomputed"] else [])
       
       
        