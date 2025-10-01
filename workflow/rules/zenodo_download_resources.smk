rule download_zenodo_resources:
    """
    Download specific files from Zenodo record.
    """
    output:
        download_complete = ZENODO_RESOURCES_DOWNLOAD_COMPLETE,
        extracted_file = directory(ZENODO_RESOURCE_DIRECTORY_UNZIPPED) # The actual file you want
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
        
        # Mark as complete
        touch {output.download_complete}
        """


rule download_zenodo_batch_analysis:
    """
    Download specific files from Zenodo record.
    """
    output:
        download_complete = ZENODO_BATCH_DOWNLOAD_COMPLETE,
        extracted_file = directory(ZENODO_BATCH_DIRECTORY_UNZIPPED) 
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
        extracted_file = directory(ZENODO_DATA_INDIV_DIRECTORY_UNZIPPED) 
    log:
        "logs/download_zenodo_pd1_itermediate_files.log"
    message:
        "Downloading PD1 intermediate files from Zenodo record ID: {params.zenodo_record_id}"
    params:
        zenodo_record_id = ZENODO_RECORD_ID,
        file = ZENODO_INDIVIDUAL_CANCERS_FILENAME  
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


rule download_combined_gdc_data:
    """
    Download specific files from Zenodo record.
    """
    output:
        download_complete = ZENODO_GDC_DOWNLOAD_COMPLETE,
        extracted_file = directory(ZENODO_GDC_DIRECTORY_UNZIPPED) 
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
        
        # Mark as complete
        touch {output.download_complete}
        """
