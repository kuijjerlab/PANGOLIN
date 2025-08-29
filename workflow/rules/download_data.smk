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