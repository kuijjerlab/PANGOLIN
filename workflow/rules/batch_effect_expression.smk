# Generate comprehensive MBatch batch effect visualization
rule plot_mbatch_results:
    """
    Create batch effect visualization from MBatch analysis results.
    """
    input:
        batch_files = BATCH_FILES
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
            --batch_files "{input.batch_files}" \
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
