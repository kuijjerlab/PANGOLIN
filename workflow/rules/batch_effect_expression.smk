
# MBatch batch effect analysis for cancer-specific expression data
rule analyze_batch_effect:
    """
    Analyze batch effects in cancer-specific normalized expression data using MBatch.
    Outputs include:
    - PCA plots showing sample clustering by batch variables
    - DSC statistics quantifying batch effect magnitude
    - Diagnostic plots for batch effect assessment
    
    Note: This analysis is computationally intensive and may require extended runtime.
    Prerequisites: MBatch conda environment (create from envs/mbatch.yaml)
    """
    input:
        expression_file = PYSNAIL_NORMALIZED_FILE_CANCER_SPECIFIC,
        group_file = GROUP_FILE,
        batch_file = BATCH_FILE,
        clinical_file = CLINICAL_FILE_RDATA 
    output:
        output_dir = directory(BATCH_DIR_CANCER)
    log:
        "logs/analyze_batch_effect_{cancer}.log"
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
            --output_directory {output.output_dir} \
        > {log} 2>&1
        """
# Generate comprehensive MBatch batch effect visualization
rule plot_mbatch_results:
    """
    Create batch effect visualization from MBatch analysis results.
    """
    input:
        batch_dir = BATCH_DIR_ALL_CANCERS
    output:
        pdf_file = BATCH_EFFECT_PDF
    log:
        "logs/plot_mbatch_results.log"
    message:
        "Creating batch effect figure"
    params:
        bin = config["bin"]
    shell:
        """
        Rscript {params.bin}/plot_mbatch_results.R \
            --batch_results_dir {input.batch_dir} \
            --output_file {output.pdf_file} \
        > {log} 2>&1
        """

# ComBat batch effect correction for multi-cancer expression data
rule correct_batch_effect:
    """
    Apply ComBat batch effect correction to normalized expression data.
    """
    input:
        expression_file = PYSNAIL_NORMALIZED_FILE,
        group_file = GROUP_FILE,
        batch_file = BATCH_FILE,
        clinical_file = CLINICAL_FILE_RDATA 
    output:
        batch_corrected_expression_file = BATCH_CORRECTED_EXPRESSION_FILE
    log:
        "logs/correct_batch_effect.log"
    message:
        "Correcting batch effect in expression data"
    params:
        bin = config["bin"]
    shell:
        """
        Rscript {params.bin}/combat_correct.R \
            --expression_file {input.expression_file} \
            --group_file {input.group_file} \
            --batch_file {input.batch_file} \
            --clin_file {input.clinical_file} \
            --output_file {output.batch_corrected_expression_file} \
        > {log} 2>&1
        """
