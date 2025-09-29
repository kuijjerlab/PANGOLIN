
# PANDA/LIONESS network inference using netZooPy
# rule run_panda_lioness:
#     """
#     Run PANDA and LIONESS network inference using netZooPy.
#     PANDA builds an initial regulatory network and LIONESS estimates 
#     sample-specific networks for each individual sample.
#     """
#     input:
#         exp_file = EXPRESSION_PANDA_FILE, 
#         motif_file = MOTIF_PANDA_FILE, 
#         ppi_file = PPI_PANDA_FILE 
#     output:
#         network_dir = directory(NETWORKS_DIR)
#     log:
#         "logs/run_panda_lioness.log"
#     params:
#         bin = config["bin"],
#         start_sample = 1,
#         end_sample = 10,  # adjust based on your sample count
#         computing = "cpu",  # change to "gpu" if you have GPU support
#         random_seed = 10,
#         ncores = 10
#     conda:
#         NETZOOPY_YAML
#     shell:
#         """
#         set +u
#         unset PYTHONPATH
#         unset PYTHONHOME
#         export PYTHONNOUSERSITE=1
#         
#         # Create networks directory
#         mkdir -p {output.network_dir}
#         
#         python {params.bin}/run_panda_lioness.py \
#             --exp_file {input.exp_file} \
#             --motif_file {input.motif_file} \
#             --ppi_file {input.ppi_file} \
#             --output_dir {output.network_dir} \
#             --start_sample {params.start_sample} \
#             --end_sample {params.end_sample} \
#             --computing {params.computing} \
#             --random_seed {params.random_seed} \
#             --ncores {params.ncores} \
#             > {log} 2>&1
#         """

# Create mapping file linking LIONESS networks to sample IDs

rule create_lioness_mapping:
    input:
        samples_file = SAMPLES_WITH_CANCER_FILE,
        networks_dir = NETWORKS_DIR
    output:
        mapping = LIONESS_SAMPLE_MAPPING
    log:
        "logs/create_lioness_sample_mapping.log"
    message:
        "Creating LIONESS sample mapping file"
    params:
        bin = config["bin"],
    shell:
        """
        Rscript {params.bin}/create_lioness_sample_mapping_file.R \
            --network_dir {input.networks_dir} \
            --samples_panda_file {input.samples_file} \
            --output_file {output.mapping} \
            > {log} 2>&1
        """


# ## Save cancer-specific LIONESS networks
rule save_cancer_specific_lioness_networks:
    """
    Save LIONESS sample-specific networks to cancer-specific RData files.
    If the number of networks is large, split them into multiple files.
    """
    input:
        network_dir = NETWORKS_DIR,
        lioness_sample_mapping = LIONESS_SAMPLE_MAPPING
    output:
        output_dir = directory(OUTPUT_DIR_FINAL_MERGED_NETWORKS)
    log:
        "logs/save_cancer_specific_lioness_networks.log"
    message:
        "Saving LIONESS networks to cancer-specific Rdata files"
    params:
        bin = config["bin"]
    shell:
        """
        Rscript {params.bin}/save_networks.R \
            --network_dir {input.network_dir} \
            --lioness_sample_mapping {input.lioness_sample_mapping} \
            --output_dir {output.output_dir} \
            > {log} 2>&1    
        """

# rule create_network_mapping:
#     """
#     Create a mapping file linking the networks RData files to its corresponding cancer type.
#     """
#     input:
#         network_merged_dir = OUTPUT_DIR_FINAL_MERGED_NETWORKS
#     output:
#         mapping_file = NETWORK_CANCER_MAPPING_FILE
#     log:
#         "logs/create_network_cancer_mapping.log"
#     message:
#         "Creating network to cancer type mapping file"
#     params:
#         bin = config["bin"]
#     shell:
#         """
#         Rscript {params.bin}/create_network_mapping.R \
#             --network_dir {input.network_merged_dir} \
#             --output_file {output.mapping_file} \
#             > {log} 2>&1
#         """

# ## Apply quantile normalization on the networks
# rule normalize_networks:
#     """
#     Apply quantile normalization to cancer-specific LIONESS networks.
#     """
#     input:
#         network_files = ALL_MERGED_NETWORKS,
#         sample_file = SAMPLES_WITH_CANCER_FILE
#     output:
#         output_dir = directory(OUTPUT_DIR_NORMALIZED_NETWORKS)
#     log:
#         "logs/quantile_normalize_networks.log"  
#     message:
#         "Applying quantile normalization to networks"
#     params:
#         bin = config["bin"]
#     shell:
#         """
#         Rscript {params.bin}/quantile_normalize_networks.R \
#             --network_files "{input.network_files}" \
#             --output_dir {output.output_dir} \
#             --sample_file {input.sample_file} \
#             > {log} 2>&1    
#         """

# ## Create a network edge file
# rule create_network_edge_file:
#     """
#     Create a network edge file from the PANDA network file.
#     """
#     input:
#         panda_input = PANDA_NETWORK_FILE
#     output:
#         edge_file = NETWORK_EDGE_FILE
#     log:
#         "logs/create_network_edge_file.log"
#     message:
#         "Creating network edge file"
#     params:
#         bin = config["bin"]
#     shell:
#         """
#         Rscript {params.bin}/create_edge_file.R \
#             --panda_network_file {input.panda_input} \
#             --output_edge_file {output.edge_file} \
#             > {log} 2>&1
#         """

# ## Calculate cancer-specific gene indegrees
# rule calculate_indegree:
#     """
#     Calculate gene indegree (number of incoming regulatory edges) from 
#     quantile-normalized LIONESS networks for each cancer type.
#     """
#     input:
#         network_dir = OUTPUT_DIR_NORMALIZED_NETWORKS,
#         edge_file = NETWORK_EDGE_FILE
#     output:
#         indegree_file = CANCER_INDEGREE_FILE
#     log:
#         "logs/calculate_indegree_{cancer}.log"
#     message:
#         "Calculating gene indegrees for cancer type: {wildcards.cancer}"
#     params:
#         bin = config["bin"]
#     shell:
#         """
#         Rscript {params.bin}/calculate_indegree.R \
#             --tumor_type {wildcards.cancer} \
#             --network_dir {input.network_dir} \
#             --edge_file {input.edge_file} \
#             --output_file {output.indegree_file} \
#             > {log} 2>&1
#         """
