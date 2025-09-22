## Extract PDL1 gene expression for each cancer type ##

rule extract_PDL1_gene_expression:
    input:
        expression_file = EXPRESSION_PANDA_FILE,
        samples_file = SAMPLES_WITH_CANCER_FILE
    output:
        out_file = OUTPUT_PDL1_EXP_CANCER
    log:
         "logs/extract_pdl1_expression_{cancer}.log"
    message:
        "Extracting expression data for PDL1: {wildcards.cancer}"
    params:
        bin = config["bin"],
        gene_id = GENE_ID # this is how it shouldbe

    shell:
        """
        Rscript {params.bin}/extract_PDL1_expression.R \
            --exp_file {input.expression_file} \
            --samples_file {input.samples_file} \
            --tumor {wildcards.cancer} \
            --gene_id {params.gene_id} \
            --output {output.out_file} \
            > {log} 2>&1
        """

## Extract PD1 pathway-based individual mappings (heterogeneity scores) ##        
rule extract_pd1_pathway_individual_mappings:
    input:
       tumor_pathways_mapping_path = TUMOR_PATHWAYS_MAPPING_PATH
    output:
        out_file = OUTPUT_CANCER_PD1_MAPPINGS
    log:
         "logs/extract_pd1_pathway_mappings_{cancer}.log"
    params:
        bin = config["bin"]
    message:
        "Extracting PD1 pathway individual mappings: {wildcards.cancer}"
    shell:
        """
        Rscript {params.bin}/extract_pd1_pathway_individual_scores.R \
            --tumor_pathways_mapping_path {input.tumor_pathways_mapping_path} \
            --output {output.out_file} \
             > {log} 2>&1
        """
