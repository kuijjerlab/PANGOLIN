libraries <- c("data.table", "dplyr", "doParallel",
            "survminer", "survival", "magrittr", "gridExtra",
            "purrr")
lapply(libraries, library, character.only = TRUE)
#' Generate bcr code from Sample ID
#'
#' This function generates a BCR code from a given sample ID.
#' The BCR code is created by concatenating the first three elements of the
#' sample ID, which are separated by hyphens.
#'
#' @param sample_id A character vector of sample ID for example 
#'                  (patient_id (f.ex TCGA-BL-A13I-01A-11R-A277-07))
#' @return A character vector containing the BCR code (f.ex TCGA-BL-A13I)
#' @export
#' 
make_bcr_code <- function(sample_id) {
        if (!is.character(sample_id)) {
            stop("sample_id must be a character vector.")
        }
        # Process sample_id and extract the first three components
        bcr_code <- sapply(strsplit(sample_id, "-"), function(split_id) {
            if (length(split_id) < 3) {
                return(NA)  # Return NA for malformed input
            }
            paste(split_id[1], split_id[2], split_id[3], sep = "-")
        })
        return(bcr_code)
    }
