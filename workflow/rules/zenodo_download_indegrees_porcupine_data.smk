rule download_zenodo_individual_cancers_directory:
    """
    Download and extract the individual cancers directory from Zenodo.
    """
    output:
        download_complete = "results/data_individual_cancers_download_complete.txt"
    log:
        "logs/download_zenodo_individual_cancers_directory.log"
    message:
        "Downloading individual cancers directory from Zenodo record ID: {params.zenodo_record_id}"
    params:
        zenodo_record_id = ZENODO_RECORD_ID,
        file = ZENODO_INDIVIDUAL_CANCERS_DIRECTORY
    shell:
        """
        # Create log directory if it doesn't exist
        mkdir -p $(dirname {log})

        # Download the tar.gz file
        echo "Downloading: {params.file}" >> {log}
        wget -O "{params.file}" "https://zenodo.org/record/{params.zenodo_record_id}/files/{params.file}" >> {log} 2>&1

        # Extract the tar.gz file
        echo "Extracting: {params.file}" >> {log}
        tar -xzvf "{params.file}" -C results/ >> {log} 2>&1

        # Mark the download as complete
        touch {output.download_complete}
        """

