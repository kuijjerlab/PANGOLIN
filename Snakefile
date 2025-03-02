## How to run this?
## snakemake --cores 1 -np ### For dry run
## module load snakemake/7.23.1-foss-2022a
## module load R-bundle-Bioconductor/3.15-foss-2022a-R-4.2.1

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

## Input Files ##
CANCER_TYPES = config["cancer_types"]
CLINICAL_FILE = config["clinical_file"]
CANCER_COLOR_FILE = config["cancer_color_file"]
TUMOR_CLIN_FILE = os.path.join(OUTPUT_DIR, "{cancer}", "clinical", "curated_clinical_{cancer}.txt")
TUMOR_PD1_DIR = os.path.join(OUTPUT_DIR, "{cancer}", "pd1_data")
TUMOR_PATHWAYS_MAPPING_PATH = os.path.join(OUTPUT_DIR, "{cancer}", "porcupine", "individual_scores_{cancer}.RData")

## Output Files ##

CANCER_LEGEND_PDF = os.path.join(FIG_DIR, "cancer_legend.pdf")
OUTPUT_CANCER = os.path.join(OUTPUT_DIR, "{cancer}", "clinical", "curated_clinical_{cancer}.txt")
OUTPUT_CANCER_COX = os.path.join(OUTPUT_DIR, "{cancer}", "cox", "{cancer}_PDL1_cox_multivariate_res.txt")
OUTPUT_CANCER_PD1_MAPPINGS  = os.path.join(OUTPUT_DIR, "{cancer}", "pd1_data", "pd1_individual_scores_norm_{cancer}.RData")
COX_RESULTS_ALL = os.path.join("data_all", "cox_results_all", "PDL1_cox_multivarite_res_all.txt")

## Parameters ##
ALPHA = config["alpha"]
NUMBER_FOLDS = config["number_folds"]
NUMBER_CORES = config["number_cores"]
NUMBER_TIMES = config["number_times"]

# Rules ##
rule all:
    input:
        expand(OUTPUT_CANCER, cancer = CANCER_TYPES),
        CANCER_LEGEND_PDF,
        expand(OUTPUT_CANCER_PD1_MAPPINGS, cancer = CANCER_TYPES),
        expand(OUTPUT_CANCER_COX, cancer = CANCER_TYPES),
        COX_RESULTS_ALL

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


rule combine_cox_results:
    input:
        expand(OUTPUT_CANCER_COX, cancer=CANCER_TYPES)
    output:
        COX_RESULTS_ALL
    message:
        "Combining all Cox regression results into one table."
    shell:
        """
        echo -e "cancer\\t$(head -n 1 {input[0]})" > {output}
        for file in {input}; do
            cancer=$(basename $(dirname $(dirname $file)))
            tail -n +2 $file | awk -v cancer=$cancer '{{print cancer"\\t"$0}}' >> {output}
        done
        """