rule download_zenodo_batch_results:
    """
    Download batch results directory from Zenodo record.
    """
    output:
        download_complete = ZENODO_BATCH_DOWNLOAD_COMPLETE
    log:
        "logs/download_zenodo_batch_results.log"
    params:
        zenodo_record_id = ZENODO_RECORD_ID,
        # Add specific filenames you want to download
        file = ZENODO_BATCH_FILENAME,
        output_dir = OUTPUT_DIR_ALL_CANCERS
    shell:
        """
        # Create directory if it doesn't exist
        mkdir -p {params.output_dir}
        
        # Download specific files
        echo "Downloading: {params.file}"
        wget -O "{params.output_dir}/{params.file}" "https://zenodo.org/records/{params.zenodo_record_id}/files/{params.file}" >> {log} 2>&1
        
        # Extract to output directory
        cd {params.output_dir}
        unzip "{params.file}" >> {log} 2>&1
        touch {output.download_complete}
        """
