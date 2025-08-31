###############################################################################
### DATA DOWNLOAD AND NORMALIZATION PIPELINE PATHS                        ###
###############################################################################

## Primary output directory for downloaded TCGA GDC data
OUTPUT_GDC_DIR = os.path.join(OUTPUT_DIR_ALL_CANCERS, "gdc_data")

## Individual cancer type expression GDC files (RData format)
## Pattern: TCGA-{cancer}.RData (e.g., TCGA-BRCA.RData, TCGA-LUAD.RData)
OUTPUT_GDC_FILE = os.path.join(OUTPUT_GDC_DIR, "TCGA-{cancer}.RData")

#------------------------------------------------------------------------------
# Combined Expression Data Processing
#------------------------------------------------------------------------------
## Directory for combined expression data from all cancers
EXPRESSION_DIR_GDC = os.path.join(OUTPUT_DIR_ALL_CANCERS, "gdc_data")

## Combined expression output directory
OUTPUT_DIR_DOWNLOAD_COMBINED = os.path.join(OUTPUT_DIR_ALL_CANCERS, "combined_gdc_data")

## Combined expression matrix file (all cancers)
## Contains log2-transformed gene expression values across all cancer types
OUTPUT_EXP_COMBINED_FILE = os.path.join(OUTPUT_DIR_DOWNLOAD_COMBINED, "hg38_STAR_counts.tsv")

## Sample group mapping file
## Maps each sample to its corresponding cancer type
GROUP_FILE = os.path.join(OUTPUT_DIR_DOWNLOAD_COMBINED, "hg38_sample_groups.tsv")

## Feature annotation file
## Contains gene annotations and genomic coordinates
FEATURE_FILE = os.path.join(OUTPUT_DIR_DOWNLOAD_COMBINED, "hg38_features.RData")

#------------------------------------------------------------------------------
# PySNAIL Normalized Expression Data
#------------------------------------------------------------------------------
## qsmooth normalized expression matrix (PySNAIL output)
## Contains batch-effect corrected log2-transformed expression values
PYSNAIL_NORMALIZED_FILE = os.path.join(OUTPUT_DIR_DOWNLOAD_COMBINED, "pysnail_normalized_STAR_counts.tsv")

#------------------------------------------------------------------------------
# Cancer-Specific Normalized Expression Files
#------------------------------------------------------------------------------
## Directory for cancer-specific normalized expression files
OUTPUT_DIR_PYSNAIL_CANCER = os.path.join(OUTPUT_DIR_ALL_CANCERS, "pysnail_normalized_individual_cancer_expression")
PYSNAIL_NORMALIZED_FILE_CANCER_SPECIFIC = os.path.join(OUTPUT_DIR_PYSNAIL_CANCER, "normalized_expression_TCGA-{cancer}.RData")

# Download TCGA GDC expression data for individual cancer types
rule download_gdc_data:
    """
    Download TCGA gene expression data from GDC for each cancer type.
    """
    output:
        out_file = OUTPUT_GDC_FILE
    message:
        "Downloading GDC expression data for cancer type: {wildcards.cancer}"
    params:
        bin = config["bin"]
    shell:
        """
        Rscript {params.bin}/download_TCGA_expression.R \
            --tumor {wildcards.cancer} \
            --output_file {output.out_file}
        """
    

# Combine TCGA expression data from all cancer types
rule combine_all_expression_data:
    """
    Combine individual cancer type expression files into one file.
    
    This rule processes downloaded TCGA expression data from multiple cancer types
    and creates the following output files:
    - Combined expression matrix (genes x samples across all cancers)
    - Sample-to-cancer group mapping file
    - Gene feature annotation file with genomic coordinates

    """
    input:
        exp_dir = EXPRESSION_DIR_GDC
    output:
        combined_expression_file = OUTPUT_EXP_COMBINED_FILE,
        group_file = GROUP_FILE,
        feature_file = FEATURE_FILE
    message:
        "Combining all expression data for all cancers"
    params:
        bin = config["bin"]
    shell:
        """
        Rscript {params.bin}/combine_allTCGA_expression.R \
            --expression_dir {input.exp_dir} \
            --combined_expression_file {output.combined_expression_file} \
            --group_file {output.group_file} \
            --feature_file {output.feature_file}
        """
# PySNAIL qsmooth normalization of combined expression data
rule qsmooth_normalization:
    """
    Apply qsmooth normalization to combined TCGA expression data using PySNAIL.
    
    Prerequisites: PySNAIL package must be installed in envs/ directory.
    Installation: git clone git@github.com:kuijjerlab/PySNAIL.git

    """
    input:
        xprs = OUTPUT_EXP_COMBINED_FILE,
        groups = GROUP_FILE
    output:
        norm = PYSNAIL_NORMALIZED_FILE
    params:
        threshold = 0.2,
        bin = config["bin"]
    conda:
        "workflow/envs/pysnail.yaml"
    shell:
        """
        python {params.bin}/normalize_with_pysnail.py {input.xprs} {input.groups} {output.norm} --threshold {params.threshold}
        """

rule split_expression_by_cancer:
    """
    Split the large normalized expression file into cancer-specific files.
    """
    input:
        expression_file = PYSNAIL_NORMALIZED_FILE,
        group_file = GROUP_FILE
    output:
        output_directory = directory(OUTPUT_DIR_PYSNAIL_CANCER)
    message:
        "Splitting expression data into cancer-specific files and saving them"
    params:
        bin = config["bin"]
    shell:
        """
        Rscript {params.bin}/save_normalized_exp_per_cancer.R \
            --expression_file {input.expression_file} \
            --group_file {input.group_file} \
            --output_dir {output.output_directory}
        """   

