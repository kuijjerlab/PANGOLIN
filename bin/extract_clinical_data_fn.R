#' Load curated clinical data from an Excel file
#'
#' @param clin_file_path Path to the clinical data Excel file.
#' @return A data frame containing the clinical data.
#' @import openxlsx
#' @export
load_clin_curated <- function(clin_file_path) {
    if (!file.exists(clin_file_path)) {
        stop("Error: The specified file does not exist: ", clin_file_path)
    }
    
    clin <- tryCatch({
        openxlsx::read.xlsx(clin_file_path, sheet = "TCGA-CDR")
    }, error = function(e) {
        stop("Error reading the Excel file: ", e$message)
    })

    # Drop first column if unnecessary (assuming it's an index)
    clin <- clin[, -1]

    return(clin)
}

#' Filter clinical data for a specific tumor type
#'
#' @param clin A data frame containing clinical data.
#' @param tumor A string representing the tumor type (e.g., "ACC").
#' @return A filtered data frame containing only the specified tumor type.
#' @export
select_tumor_clin_curated <- function(clin, tumor) {
    if (!is.character(tumor) || length(tumor) != 1 || tumor == "") {
        stop("Error: The 'tumor' argument must be a non-empty string.")
    }

    clin_tumor <- clin[grepl(tumor, clin$type, ignore.case = TRUE), ]

    if (nrow(clin_tumor) == 0) {
        warning("No matches found for tumor type: ", tumor)
    }

    return(clin_tumor)
}


