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

# Always include all rule files (Snakemake will only execute needed rules)
include: "workflow/rules/zenodo_download_resources.smk"
include: "workflow/rules/batch_effect_expression.smk"
include: "workflow/rules/download_normalize_expression_data.smk"
include: "workflow/rules/zenodo_download_batch_results.smk"


# Rules ##
rule all:
    input:
        # Always download resources (both workflows need this)
        ZENODO_RESOURCES_DOWNLOAD_COMPLETE,
        
        # Conditional inputs based on analysis type

        ([ZENODO_BATCH_DOWNLOAD_COMPLETE] if ANALYSIS_TYPE == "precomputed" else []) +
        
        # Full workflow specific outputs
        ([expand(OUTPUT_GDC_FILE, cancer=CANCER_TYPES)] if ANALYSIS_TYPE == "full_workflow" else [])
        # ([OUTPUT_EXP_COMBINED_FILE] if ANALYSIS_TYPE == "full_workflow" else []) 
        # ([GROUP_FILE] if ANALYSIS_TYPE == "full_workflow" else []) +
        # ([FEATURE_FILE] if ANALYSIS_TYPE == "full_workflow" else []) +
        # ([PYSNAIL_NORMALIZED_FILE] if ANALYSIS_TYPE == "full_workflow" else []) +
        # ([OUTPUT_DIR_PYSNAIL_CANCER] if ANALYSIS_TYPE == "full_workflow" else []) +
        # ([expand(BATCH_DIR_CANCER, cancer=CANCER_TYPES)] if ANALYSIS_TYPE == "full_workflow" else []) +
        # ([BATCH_CORRECTED_EXPRESSION_FILE] if ANALYSIS_TYPE == "full_workflow" else []) +
        
        # # Shared outputs (both workflows can generate this)
        # ([BATCH_EFFECT_PDF] if ANALYSIS_TYPE in ["full_workflow", "precomputed"] else [])


