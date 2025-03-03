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
PPI_FILE = config["ppi_file"]
MOTIF_FILE = config["motif_file"]

## Output Files ##

CANCER_LEGEND_PDF = os.path.join(FIG_DIR, "cancer_legend.pdf")
OUTPUT_CANCER = os.path.join(OUTPUT_DIR, "{cancer}", "clinical", "curated_clinical_{cancer}.txt")

## Output Files for Univariate Cox on PD1-pathway based heterogeneity scores ##
OUTPUT_CANCER_UNIVARIATE_COX_SUMMARY = os.path.join(OUTPUT_DIR, "{cancer}", "cox", "{cancer}_PD1_pathway_cox_univariate_model_summary.txt")
OUTPUT_CANCER_UNIVARIATE_COX_PREDICTED_SCORES = os.path.join(OUTPUT_DIR, "{cancer}", "cox", "{cancer}_PD1_pathway_cox_univariate_predited_risk_scores.txt")
UNIVARIATE_COX_SUMMARY_ALL = os.path.join("data_all", "cox_results_all", "PD1_pathway_cox_univariate_model_summary_all.txt")
UNIVARIATE_COX_PREDICTED_SCORES_ALL = os.path.join("data_all", "cox_results_all", "PD1_pathway_cox_univariate_predited_risk_scores_all.txt")


## Output Files for multivariate regularized Cox on PDL1-edges ##
OUTPUT_CANCER_PD1_MAPPINGS  = os.path.join(OUTPUT_DIR, "{cancer}", "pd1_data", "pd1_individual_scores_norm_{cancer}.RData")
OUTPUT_CANCER_COX = os.path.join(OUTPUT_DIR, "{cancer}", "cox", "{cancer}_PDL1_cox_multivariate_res.txt")
COX_RESULTS_ALL = os.path.join("data_all", "cox_results_all", "PDL1_cox_multivarite_res_all.txt")
PDL1_CIRCULAR_PLOT = os.path.join(FIG_DIR, "circular_pdl1_plot_{threshold_cox}.pdf")


## Parameters ##
ALPHA = config["alpha"]
NUMBER_FOLDS = config["number_folds"]
NUMBER_CORES = config["number_cores"]
NUMBER_TIMES = config["number_times"]
THESHOLD_COX = config["threshold_cox"]

# Rules ##
rule all:
    input:
        expand(OUTPUT_CANCER, cancer = CANCER_TYPES),
        CANCER_LEGEND_PDF,
        expand(OUTPUT_CANCER_PD1_MAPPINGS, cancer = CANCER_TYPES),
        expand(OUTPUT_CANCER_UNIVARIATE_COX_SUMMARY, cancer = CANCER_TYPES),
        expand(OUTPUT_CANCER_UNIVARIATE_COX_PREDICTED_SCORES, cancer = CANCER_TYPES),
        UNIVARIATE_COX_SUMMARY_ALL,
        UNIVARIATE_COX_PREDICTED_SCORES_ALL,
        expand(OUTPUT_CANCER_COX, cancer = CANCER_TYPES),
        COX_RESULTS_ALL, 
        expand(PDL1_CIRCULAR_PLOT, threshold_cox = THESHOLD_COX)

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
## Extract PD1 pathway individual mappings ##        
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
        bin = config["bin"],
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
        COX_RESULTS_ALL
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
        cox_results_all = COX_RESULTS_ALL,
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
            --cox_results_all {input.cox_results_all} \
            --ppi_file {input.ppi_file} \
            --motif_file {input.motif_file} \
            --threshold {wildcards.threshold_cox} \
            --output {output.out_file}
        """