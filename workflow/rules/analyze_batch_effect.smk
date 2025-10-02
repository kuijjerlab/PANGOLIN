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
        batch_files = expand(
            os.path.join(BATCH_DIR_ALL_CANCERS, "TCGA-{{cancer}}", "{batch}", "ManyToMany/1.0/1.0/PCAAnnotations.tsv"),
            batch=BATCHES
        )
    log:
        "logs/analyze_batch_effect_{cancer}.log"
    message:
        "Analyzing batch effect for: {wildcards.cancer}"
    params:
        bin = config["bin"],
        nthreads = config["number_cores_mbatch"],
        output_dir = BATCH_DIR_CANCER
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
            --output_directory {params.output_dir} \
        > {log} 2>&1
        """