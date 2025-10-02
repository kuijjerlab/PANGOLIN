
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
        
