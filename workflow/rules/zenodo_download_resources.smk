# Rule ordering to resolve ambiguity between download and compute rules
# For precomputed analysis, download rules should take precedence
# For full workflow, compute rules should take precedence
ruleorder: combine_all_expression_data > download_combined_gdc_data

rule download_zenodo_resources:
    """
    Download specific files from Zenodo record.
    """
    priority: 100
    output:
        download_complete = ZENODO_RESOURCES_DOWNLOAD_COMPLETE,
        extracted_file = directory(ZENODO_RESOURCE_DIRECTORY_UNZIPPED),
        clinical_file = CLINICAL_FILE,
        cancer_color_file = CANCER_COLOR_FILE,
        pathways_hierarchy = PATHWAYS_HIERARCHY_FILE,
        pathways_hsa_id_file = PATHWAYS_HSA_ID_FILE,
        list_pathways = LIST_PATHWAYS_FILE,
        gmt_file = GMT_FILE,
        ppi_file = PPI_FILE,
        motif_file = MOTIF_FILE,
        samples_file = SAMPLES_FILE,
        immune_file = IMMUNE_FILE,
        batch_file = BATCH_FILE,
        clinical_rdata = CLINICAL_FILE_RDATA,
        summary_table_pd1 = SUMMARY_TABLE_PD1

    log:
        "logs/download_zenodo_resources.log"
    message:
        "Downloading resources from Zenodo record ID: {params.zenodo_record_id}"
    params:
        zenodo_record_id = ZENODO_RECORD_ID,
        file = ZENODO_RESOURCE_FILENAME  
    shell:
        """
        # Download the zip file
        echo "Downloading: {params.file}"
        wget -O "{params.file}" "https://zenodo.org/records/{params.zenodo_record_id}/files/{params.file}?download=1" >> {log} 2>&1
        
        # Unzip the file, excluding __MACOSX folders
        echo "Extracting: {params.file}"
        unzip -o "{params.file}" -x "__MACOSX/*" >> {log} 2>&1
        
        # Remove any __MACOSX folders that might have been created anyway
        find . -name "__MACOSX" -type d -exec rm -rf {{}} + 2>/dev/null || true
        
        # Remove the zip file after extraction
        echo "Removing zip file: {params.file}"
        rm -f "{params.file}" >> {log} 2>&1
        
        # Verify all expected files exist
        echo "Verifying extracted files..."
        for file in {output.clinical_file} {output.cancer_color_file} {output.summary_table_pd1} {output.clinical_rdata} {output.batch_file} {output.pathways_hierarchy} {output.pathways_hsa_id_file} {output.list_pathways} {output.gmt_file} {output.ppi_file} {output.motif_file} {output.samples_file} {output.immune_file}; do
            if [ ! -f "$file" ]; then
                echo "Error: Expected file $file not found after extraction" >> {log}
                exit 1
            fi
        done
        echo "All expected files verified successfully"
        
        # Mark as complete
        touch {output.download_complete}
        """


rule download_zenodo_batch_analysis:
    """
    Download specific files from Zenodo record.
    """
    priority: 100
    output:
        download_complete = ZENODO_BATCH_DOWNLOAD_COMPLETE,
        extracted_file = directory(ZENODO_BATCH_DIRECTORY_UNZIPPED),
        batch_files = BATCH_FILES
    log:
        "logs/download_zenodo_batch_analysis.log"
    message:
        "Downloading batch_analysis from Zenodo record ID: {params.zenodo_record_id}"
    params:
        zenodo_record_id = ZENODO_RECORD_ID,
        file = ZENODO_BATCH_FILENAME  
    shell:
        """
        # Create target directory if it doesn't exist
        mkdir -p results/data_all >> {log} 2>&1
        
        # Download the zip file
        echo "Downloading: {params.file}"
        wget -O "{params.file}" "https://zenodo.org/records/{params.zenodo_record_id}/files/{params.file}?download=1" >> {log} 2>&1
        
        # Unzip the file directly to target directory, excluding __MACOSX folders
        echo "Extracting: {params.file} to results/data_all/"
        unzip -o "{params.file}" -d "results/data_all/" -x "__MACOSX/*" >> {log} 2>&1
        
        # Remove any __MACOSX folders that might have been created anyway
        find results/data_all/ -name "__MACOSX" -type d -exec rm -rf {{}} + 2>/dev/null || true
        
        # Ensure the expected output directory exists (Snakemake expects this)
        mkdir -p "{output.extracted_file}" >> {log} 2>&1
        
        # Remove the zip file after extraction
        echo "Removing zip file: {params.file}"
        rm -f "{params.file}" >> {log} 2>&1
        
        # Mark as complete
        touch {output.download_complete}
        """


rule download_zenodo_pd1_itermediate_files:
    """
    Download specific files from Zenodo record.
    """
    output:
        download_complete = ZENODO_PD1_DOWNLOAD_COMPLETE,
        extracted_file = directory(ZENODO_DATA_INDIV_DIRECTORY_UNZIPPED),
        pd1_links = expand(TUMOR_PD1_LINKS, cancer=CANCER_TYPES),
        pd1_net = expand(TUMOR_PD1_NET, cancer=CANCER_TYPES),
        indegree = expand(CANCER_INDEGREE_FILE, cancer = CANCER_TYPES)
    log:
        "logs/download_zenodo_pd1_itermediate_files.log"
    message:
        "Downloading PD1 intermediate files from Zenodo record ID: {params.zenodo_record_id}"
    params:
        zenodo_record_id = ZENODO_RECORD_ID,
        file = ZENODO_INDIVIDUAL_CANCERS_FILENAME,
 
    shell:
        """
        # Create target directory if it doesn't exist
        mkdir -p results/data_all >> {log} 2>&1
        
        # Download the zip file
        echo "Downloading: {params.file}"
        wget -O "{params.file}" "https://zenodo.org/records/{params.zenodo_record_id}/files/{params.file}?download=1" >> {log} 2>&1
        
        # Unzip the file directly to target directory, excluding __MACOSX folders
        echo "Extracting: {params.file} to results/data_all/"
        unzip -o "{params.file}" -d "results/data_all/" -x "__MACOSX/*" >> {log} 2>&1
        
        # Remove any __MACOSX folders that might have been created anyway
        find results/data_all/ -name "__MACOSX" -type d -exec rm -rf {{}} + 2>/dev/null || true
        
        # Ensure the expected output directory exists (Snakemake expects this)
        mkdir -p "{output.extracted_file}" >> {log} 2>&1
        
        # Remove the zip file after extraction
        echo "Removing zip file: {params.file}"
        rm -f "{params.file}" >> {log} 2>&1
        
        # Mark as complete
        touch {output.download_complete}
        """

# Rule ordering to resolve ambiguity between download and compute rules
# For precomputed analysis, download rules should take precedence
# For full workflow, compute rules should take precedence
ruleorder: combine_all_expression_data > download_combined_gdc_data

rule download_combined_gdc_data:
    """
    Download specific files from Zenodo record.
    """
    output:
        download_complete = ZENODO_GDC_DOWNLOAD_COMPLETE,
        extracted_file = directory(ZENODO_GDC_DIRECTORY_UNZIPPED),
        output_exp_combined_file = OUTPUT_EXP_COMBINED_FILE,
        group_file = GROUP_FILE,
        feature_file = FEATURE_FILE 
    log:
        "logs/download_zenodo_pd1_itermediate_files.log"
    message:
        "Downloading GDC files from Zenodo record ID: {params.zenodo_record_id}"
    params:
        zenodo_record_id = ZENODO_RECORD_ID,
        file = ZENODO_GDC_FILENAME  
    shell:
        """
        # Create target directory if it doesn't exist
        mkdir -p results/data_all >> {log} 2>&1
        
        # Download the zip file
        echo "Downloading: {params.file}"
        wget -O "{params.file}" "https://zenodo.org/records/{params.zenodo_record_id}/files/{params.file}?download=1" >> {log} 2>&1
        
        # Unzip the file directly to target directory, excluding __MACOSX folders
        echo "Extracting: {params.file} to results/data_all/"
        unzip -o "{params.file}" -d "results/data_all/" -x "__MACOSX/*" >> {log} 2>&1
        
        # Remove any __MACOSX folders that might have been created anyway
        find results/data_all/ -name "__MACOSX" -type d -exec rm -rf {{}} + 2>/dev/null || true
        
        # Ensure the expected output directory exists (Snakemake expects this)
        mkdir -p "{output.extracted_file}" >> {log} 2>&1
        
        # Remove the zip file after extraction
        echo "Removing zip file: {params.file}"
        rm -f "{params.file}" >> {log} 2>&1
        
        # Verify all expected files exist
        echo "Verifying extracted files..."
        for file in {output.output_exp_combined_file} {output.group_file} {output.feature_file}; do
            if [ ! -f "$file" ]; then
                echo "Error: Expected file $file not found after extraction" >> {log}
                exit 1
            fi
        done
        echo "All expected files verified successfully"

        # Mark as complete
        touch {output.download_complete}
        """
