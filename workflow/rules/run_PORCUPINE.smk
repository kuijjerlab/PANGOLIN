rule run_PORCUPINE:
    input:
        network_file = CANCER_NETWORK_NORMALIZED_FILE,
        edge_file = NETWORK_EDGE_FILE,
        pathway_file = GMT_FILE
    output:
        pathways_results = PORCUPINE_PATHWAYS_RESULTS,
        pathways_results_random = PORCUPINE_PATHWAYS_RESULTS_RANDOM,
        porcupine_results = PORCUPINE_RESULTS,
        individual_scores = INDIVIDUAL_SCORES
    log:
        "logs/run_PORCUPINE_{cancer}.log"
    message:
        "Running PORCUPINE for {wildcards.cancer}"
    params:
        bin = config["bin"],
        ncores_porcupine = NCORES_PORCUPINE
    shell:
        """
        Rscript {params.bin}/runPORCUPINE_norm.R \
            --tumor_type  {wildcards.cancer} \
            --network_file {input.network_file} \
            --edge_file {input.edge_file} \
            --pathway_file {input.pathway_file} \
            --ncores {params.ncores_porcupine} \
            --pathways_results {output.pathways_results} \
            --pathways_results_random {output.pathways_results_random} \
            --porcupine_results {output.porcupine_results} \
            --individual_scores {output.individual_scores} \
            > {log} 2>&1
        """

rule extract_PD1_information:
    input:
        network_files = CANCER_NETWORK_NORMALIZED_FILE,
        edge_file = NETWORK_EDGE_FILE,
        pathway_file = GMT_FILE
    output:
        pd1_edges = TUMOR_PD1_LINKS,
        pd1_net = TUMOR_PD1_NET,
    log:
        "logs/extract_PD1_information_{cancer}.log"
    message:
        "Extracting PD-1 related information for {wildcards.cancer}"
    params:
        bin = config["bin"]
    shell:
        """
        Rscript {params.bin}/extract_pd1_information.R \
            --network_file "{input.network_files}" \
            --edge_file {input.edge_file} \
            --pathway_file {input.pathway_file} \
            --tumor_type {wildcards.cancer} \
            --pd1_edges_file {output.pd1_edges} \
            --pd1_net_file {output.pd1_net} \
            > {log} 2>&1
        """