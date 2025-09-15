## PREPARE FOR PANDA
rule prepare_files_for_PANDA:
    """
    Prepare expression data, motif priors, and PPI priors for PANDA analysis.
    Filters for protein coding genes, removes duplicates, applies expression
    thresholds, and creates PANDA-compatible input files.
    """
    input:
        expression_file = BATCH_CORRECTED_EXPRESSION_FILE,
        batch_file = BATCH_FILE,
        motif_file = MOTIF_FILE,
        ppi_file = PPI_FILE,
        samples_file = SAMPLES_FILE,
        feature_file = FEATURE_FILE,
        groups_file = GROUP_FILE,
    output:
        motif_filtered = MOTIF_PANDA_FILE,
        ppi_filtered = PPI_PANDA_FILE,
        expression_filtered = EXPRESSION_PANDA_FILE,
        samples_filtered = SAMPLES_PANDA_FILE,
        samples_with_cancer = SAMPLES_WITH_CANCER_FILE
    log:
        "logs/prepare_for_PANDA.log"
    params:
        bin = config["bin"],
        min_sample_expression = 20
    shell:
        """
        Rscript {params.bin}/prepare_for_PANDA.R \
            --expression_file {input.expression_file} \
            --batch_file {input.batch_file} \
            --motif_file {input.motif_file} \
            --ppi_file {input.ppi_file} \
            --samples_file {input.samples_file} \
            --feature_file {input.feature_file} \
            --groups_file {input.groups_file} \
            --output_motif_file_filtered {output.motif_filtered} \
            --output_ppi_file_filtered {output.ppi_filtered} \
            --output_expression_file_filtered {output.expression_filtered} \
            --output_samples_file_filtered {output.samples_filtered} \
            --output_samples_file_filtered_with_cancer_type {output.samples_with_cancer} \
            --min_sample_expression {params.min_sample_expression} \
            > {log} 2>&1
        """


# PANDA/LIONESS network inference using netZooPy
rule run_panda_lioness:
    """
    Run PANDA and LIONESS network inference using netZooPy.
    PANDA builds an initial regulatory network and LIONESS estimates 
    sample-specific networks for each individual sample.
    """
    input:
        exp_file = EXPRESSION_PANDA_FILE, 
        motif_file = MOTIF_PANDA_FILE, 
        ppi_file = PPI_PANDA_FILE 
    output:
        network_dir = directory(NETWORKS_DIR)
    log:
        "logs/run_panda_lioness.log"
    params:
        bin = config["bin"],
        start_sample = 1,
        end_sample = 10,  # adjust based on your sample count
        computing = "cpu",  # change to "gpu" if you have GPU support
        random_seed = 10,
        ncores = 10
    conda:
        "workflow/envs/netzoopy-local.yaml"
    shell:
        """
        set +u
        unset PYTHONPATH
        unset PYTHONHOME
        export PYTHONNOUSERSITE=1
        
        # Create networks directory
        mkdir -p {output.network_dir}
        
        python {params.bin}/run_panda_lioness.py \
            --exp_file {input.exp_file} \
            --motif_file {input.motif_file} \
            --ppi_file {input.ppi_file} \
            --output_dir {output.network_dir} \
            --start_sample {params.start_sample} \
            --end_sample {params.end_sample} \
            --computing {params.computing} \
            --random_seed {params.random_seed} \
            --ncores {params.ncores} \
            > {log} 2>&1
        """

