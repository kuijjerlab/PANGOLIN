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
        NETZOOPY_YAML
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

## Create LIONESS sample-to-network mapping file
rule create_lioness_mapping_file:
    """
    Create a comprehensive mapping file that associates each LIONESS
    sample-specific network file with its corresponding sample ID and cancer type.
    """
    input:
        network_dir = NETWORKS_DIR,
        samples_file = SAMPLES_WITH_CANCER_FILE
    output:
        mapping_file = LIONESS_SAMPLE_MAPPING
    log:
        "logs/create_lioness_sample_mapping_file.log"
    message:
        "Creating LIONESS sample-to-network mapping file from {input.network_dir}"
    params:
        bin = config["bin"]
    shell:
        """
        Rscript {params.bin}/create_lioness_sample_mapping_file.R \
            --network_dir {input.network_dir} \
            --samples_panda_file {input.samples_file} \
            --output_file {output.mapping_file} \
            > {log} 2>&1
        """


## Save cancer-specific LIONESS networks
rule save_cancer_specific_lioness_networks:
    """
    Save LIONESS sample-specific networks to cancer-specific RData files.
    If the number of networks is large, split them into multiple files.
    """
    input:
        network_dir = NETWORKS_DIR,
        lioness_sample_mapping = LIONESS_SAMPLE_MAPPING
    output:
        output_dir_merged = directory(OUTPUT_DIR_FINAL_MERGED_NETWORKS)
    log:
        "logs/save_cancer_specific_lioness_networks.log"
    message:
        "Saving cancer-specific LIONESS networks"
    params:
        bin = config["bin"]
    shell:
        """
         Rscript {params.bin}/save_networks.R \
            --network_dir {input.network_dir} \
            --lioness_sample_mapping {input.lioness_sample_mapping} \
            --output_dir {output.output_dir_merged} \
            > {log} 2>&1    
        """


## Apply quantile normalization on the networks
rule normalize_networks:
    """
    Apply quantile normalization to cancer-specific LIONESS networks.
    """
    input:
        network_dir = OUTPUT_DIR_FINAL_MERGED_NETWORKS,
        sample_file = SAMPLES_WITH_CANCER_FILE
    output:
        output_dir = directory(OUTPUT_DIR_NORMALIZED_NETWORKS)
    log:
        "logs/quantile_normalize_networks.log"  
    message:
        "Applying quantile normalization to networks"
    params:
        bin = config["bin"]
    shell:
        """
        Rscript {params.bin}/quantile_normalize_networks.R \
            --network_dir {input.network_dir} \
            --output_dir {output.output_dir} \
            --sample_file {input.sample_file} \
            > {log} 2>&1    
        """

## Create a network edge file
rule create_network_edge_file:
    """
    Create a network edge file from the PANDA network file.
    """
    input:
        panda_input = PANDA_NETWORK_FILE
    output:
        edge_file = NETWORK_EDGE_FILE
    log:
        "logs/create_network_edge_file.log"
    message:
        "Creating network edge file"
    params:
        bin = config["bin"]
    shell:
        """
        Rscript {params.bin}/create_edge_file.R \
            --panda_network_file {input.panda_input} \
            --output_edge_file {output.edge_file} \
            > {log} 2>&1
        """

## Calculate cancer-specific gene indegrees
rule calculate_indegree:
    """
    Calculate gene indegree (number of incoming regulatory edges) from 
    quantile-normalized LIONESS networks for each cancer type.
    """
    input:
        network_dir = OUTPUT_DIR_NORMALIZED_NETWORKS,
        edge_file = NETWORK_EDGE_FILE
    output:
        indegree_file = CANCER_INDEGREE_FILE
    log:
        "logs/calculate_indegree_{cancer}.log"
    message:
        "Calculating gene indegrees for cancer type: {wildcards.cancer}"
    params:
        bin = config["bin"]
    shell:
        """
        Rscript {params.bin}/calculate_indegree.R \
            --tumor_type {wildcards.cancer} \
            --network_dir {input.network_dir} \
            --edge_file {input.edge_file} \
            --output_file {output.indegree_file} \
            > {log} 2>&1
        """
