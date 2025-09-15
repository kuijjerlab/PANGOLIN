rule download_zenodo_resources:
    """
    Download specific files from Zenodo record.
    """
    output:
        download_complete = ZENODO_RESOURCES_DOWNLOAD_COMPLETE
    log:
        "logs/download_zenodo_resources.log"
    message:
        "Downloading resources from Zenodo record ID: {params.zenodo_record_id}"
    params:
        zenodo_record_id = ZENODO_RECORD_ID,
        # Add specific filenames you want to download
        file = ZENODO_RESOURCE_FILENAME  
    shell:
        """
        # Download specific files
        echo "Downloading: {params.file}"
        wget -O "{params.file}" "https://zenodo.org/records/{params.zenodo_record_id}/files/{params.file}" >> {log} 2>&1
        unzip "{params.file}" >> {log} 2>&1
        touch {output.download_complete}
        """

