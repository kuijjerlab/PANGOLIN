
# Download TCGA GDC expression data for individual cancer types
rule download_gdc_data:
    """
    Download TCGA gene expression data from GDC for each cancer type.
    """
    output:
        out_file = OUTPUT_GDC_FILE
    log:
        "logs/download_gdc_expression_{cancer}.log"
    message:
        "Downloading GDC expression data for cancer type: {wildcards.cancer}"
    params:
        bin = config["bin"]
    shell:
        """
        Rscript {params.bin}/download_TCGA_expression.R \
            --tumor {wildcards.cancer} \
            --output_file {output.out_file} \
        > {log} 2>&1
        """
    


rule combine_all_expression_data:
    input:
        gdc_files = expand(OUTPUT_GDC_FILE, cancer=CANCER_TYPES)
    output:
        combined_expression = OUTPUT_EXP_COMBINED_FILE,
        sample_groups = GROUP_FILE,
        features = FEATURE_FILE
    log:
        "logs/combine_all_TCGA_expression.log"
    message:
        "Combining all TCGA cancer expression data into a single matrix"
    shell:
        """
        Rscript workflow/bin/combine_allTCGA_expression.R \
            --gdc_files "{input.gdc_files}" \
            --combined_expression_file {output.combined_expression} \
            --group_file {output.sample_groups} \
            --feature_file {output.features} \
        > {log} 2>&1
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
    log:
        "logs/qsmooth_normalization.log"
    message:
        "Performing qsmooth normalization on combined expression data"
    params:
        threshold = 0.2,
        bin = config["bin"]
    conda:
        PYSNAIL_YAML_RELATIVE 
    shell:
        """
        set +u
        unset PYTHONPATH
        unset PYTHONHOME
        export PYTHONNOUSERSITE=1
        python {params.bin}/normalize_with_pysnail.py {input.xprs} {input.groups} {output.norm} --threshold {params.threshold} \
        > {log} 2>&1
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
    log:
        "logs/split_expression_by_cancer.log"
    message:
        "Splitting expression data into cancer-specific files and saving them"
    params:
        bin = config["bin"]
    shell:
        """
        Rscript {params.bin}/save_normalized_exp_per_cancer.R \
            --expression_file {input.expression_file} \
            --group_file {input.group_file} \
            --output_dir {output.output_directory} \
        > {log} 2>&1
        """

