
## Zenodo configuration ##
ZENODO_RECORD_ID = config["zenodo_record_id"]
ZENODO_RESOURCE_FILENAME = config["zenodo_resource_filename"]

##############################################################################
### ZENODO RESOURCE DOWNLOAD PATHS                                        ###
###############################################################################

## Marker file to indicate successful download of Zenodo resources
ZENODO_DOWNLOAD_COMPLETE = ".zenodo_download_complete"


rule download_zenodo_resources:
    """
    Download specific files from Zenodo record.
    """
    output:
        download_complete = ZENODO_DOWNLOAD_COMPLETE
    params:
        zenodo_record_id = ZENODO_RECORD_ID,
        # Add specific filenames you want to download
        file = ZENODO_RESOURCE_FILENAME  
    shell:
        """
        # Download specific files
        echo "Downloading: {params.file}"
        wget -O "{params.file}" "https://zenodo.org/records/{params.zenodo_record_id}/files/{params.file}"
        unzip "{params.file}"
        touch {output.download_complete}
        """

