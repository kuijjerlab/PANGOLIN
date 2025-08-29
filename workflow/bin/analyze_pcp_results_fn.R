#' Filter Porcupine Results
#'
#' This function filters the results of Porcupine analysis based on specified 
#' thresholds for adjusted p-values, effect size, and variance.
#'
#' @param pcp_res_file A string representing the file path to the PCP results.
#' @param var_th A numeric value specifying the variance threshold 
#'   (default is 25).
#' @param padjust_th A numeric value specifying the adjusted p-value threshold 
#'   (default is 0.05).
#' @param es_quant A numeric value representing the quantile cutoff for the 
#'   effect size (default is 0.95).
#' @return A data frame containing the filtered and sorted results.
#' @import data.table
#' @export
#' 
filter_pcp_results <- function(pcp_res_file,
        var_th = 25,
        padjust_th = 0.05,
        es_quant = 0.95) {

        if (!file.exists(pcp_res_file)) {
                stop("Error: File does not exist - ", pcp_res_file)
        }
        res <- tryCatch({
                fread(pcp_res_file)
        }, error = function(e) {
                stop("Error reading the file: ", e$message)
        })
        required_cols <- c("p.adjust", "es", "pc1")
        if (!all(required_cols %in% names(res))) {
                missing_cols <- setdiff(required_cols, names(res))
                stop("Error: Missing required columns: ", 
                                    paste(missing_cols, collapse = ", "))
        }

        res <- res[!is.na(res$p.adjust) & !is.na(res$es) & !is.na(res$pc1), ]
        res <- res[res$p.adjust <= padjust_th, ]
        es_cutoff <- quantile(res$es, es_quant, na.rm = TRUE)
        res <- res[res$es >= es_cutoff & res$pc1 >= var_th, ]
        res <- res[order(res$es, decreasing = TRUE), ]
        return(res)
}

#' Intersect Pathways
#'
#' This function calculates the intersection of pathways from a list of pathway
#' result sets. The intersection can be calculated either by the number of 
#' common elements or by the Jaccard index.
#'
#' @param res_list A list of pathway result sets, where each element is a 
#' vector of pathways.
#' @param type A character string specifying the type of intersection 
#' calculation. It can be either "number" for the number of common pathways 
#' or "jaccard" for the Jaccard index. Default is "number".
#' @return A matrix containing the intersection values between each pair of 
#' pathway result sets.
#' @export
intersect_pathways <- function(res_list, type = c("number", "jaccard")) {
        if (!is.list(res_list) || length(res_list) == 0) {
                stop("Error: 'res_list' must be a non-empty list of pathways.")
        }

        type <- match.arg(type)

        n <- length(res_list)
        intersection_matrix <- matrix(0, nrow = n, ncol = n)
        for (i in 1:n) {
                for (j in 1:n) {
                        if (type == "number") {
                                intersection_matrix[i, j] <- length(
                                        intersect(res_list[[i]], res_list[[j]])
                                )
                        } else if (type == "jaccard") {
                                intersection_matrix[i, j] <- jaccard(
                                        res_list[[i]], res_list[[j]]
                                )
                        }
                }
        }

        rownames(intersection_matrix) <- names(res_list)
        colnames(intersection_matrix) <- names(res_list)

        return(intersection_matrix)
}
#' Calculate the Jaccard Index
#'
#' This function calculates the Jaccard index, a measure of similarity between
#' two sets. The Jaccard index is defined as the size of the intersection 
#' divided by the size of the union of the sets.
#'
#' @param a A vector representing the first set.
#' @param b A vector representing the second set.
#' @return A numeric value representing the Jaccard index, which ranges from 
#' 0 (no similarity) to 1 (complete similarity).
#' @export
jaccard <- function(a, b) {
        if (!is.vector(a) || !is.vector(b)) {
                stop("Error: Both inputs must be vectors.")
        }

        if (length(a) == 0 || length(b) == 0) {
                stop("Error: Both input sets must be non-empty.")
        }

        intersection <- length(intersect(a, b))
        union <- length(unique(c(a, b)))
        return(intersection / union)
}


# #' Add variance (PC1) information to PCP results
# #'
# #' This function reads PCP and PCA results from specified file paths, matches 
# #' pathways between the results, and adds PC1 variance information from the 
# #' PCA results to the PCP results.
# #'
# #' @param pcp_file A character string specifying the file path to the PCP 
# #'                 results file.
# #' @param pca_file A character string specifying the file path to the PCA 
# #'                 results file (pathways_results).
# #'
# #' @return A data.table of PCP results with an added "pc1" column containing 
# #'         variance information from the PCA results.
# #' @import data.table
# #' @export


# add_variance_to_pcp_results <- function(pcp_file, pca_file) {
#         # Validate inputs
#         if (!is.character(pcp_file) || length(pcp_file) != 1) {
#                 stop("'pcp_file' must be a single character string.")
#         }
        
#         if (!is.character(pca_file) || length(pca_file) != 1) {
#                 stop("'pca_file' must be a single character string.")
#         }

#         if (!file.exists(pcp_file)) {
#                 stop(paste("PCP file does not exist:", pcp_file))
#         }
        
#         if (!file.exists(pca_file)) {
#                 stop(paste("PCA file does not exist:", pca_file))
#         }

#         # Read in PCP and PCA results
#         pcp_res <- tryCatch({
#                 fread(pcp_file)
#         }, error = function(e) {
#                 stop(paste("Error reading PCP file:", e$message))
#         })

#         pca_res <- tryCatch({
#                 fread(pca_file)
#         }, error = function(e) {
#                 stop(paste("Error reading PCA file:", e$message))
#         })

#         # Validate required columns exist
#         if (!"pathway" %in% names(pcp_res)) {
#                 stop("PCP results file must contain a 'pathway' column.")
#         }
        
#         if (!"pathway" %in% names(pca_res)) {
#                 stop("PCA results file must contain a 'pathway' column.")
#         }
        
#         if (!"pc1" %in% names(pca_res)) {
#                 stop("PCA results file must contain a 'pc1' column.")
#         }

#         # Add PC1 variance information to PCP results
#         pcp_res$pc1 <- pca_res$pc1[match(pcp_res$pathway, pca_res$pathway)]

#         # Check for unmatched pathways
#         if (any(is.na(pcp_res$pc1))) {
#                 warning("Some pathways in PCP results were not matched to PCA.")
#         }

#         return(pcp_res)
# }


#' Assign Functional Categories to Pathways
#'
#' This function assigns functional categories to pathways by integrating data 
#' from MSigDB-based pathway ID, hierarchy, and title files. It matches 
#' selected pathways to their respective hierarchy levels and provides 
#' summarized category assignments.
#'
#' @param pathways_hsa_ids_file A string representing the path to the file 
#' containing mappings of pathway names to their corresponding HSA IDs.
#' @param pathways_hierachy_file A string representing the path to the file 
#' defining the hierarchical relationships between pathways (e.g., parent-child).
#' @param list_of_pathways_file A string representing the path to the file 
#' containing a list of pathways with organism annotations.
#' @param pathways A data.table containing selected pathways with a `pathway` 
#' column to be mapped and categorized.
#'
#' @return A data.table summarizing the input pathways along with their 
#' assigned parent categories and hierarchy levels.
#'
#' @import data.table
#' @export
#' 
assign_functions_to_pathways <- function(pathways_hsa_ids_file,
                                        pathways_hierarchy_file,
                                        list_of_pathways_file,
                                        pathways) {
        ptw_ids <- load_pathway_ids(pathways_hsa_ids_file)
        ptw_hr <- load_pathway_hierarchy(pathways_hierarchy_file)
        ptw_titles <- load_pathway_titles(list_of_pathways_file)
        ptw_titles <- ptw_titles[ptw_titles$V3 == "Homo sapiens", ]
        colnames(ptw_hr) <- c("parent_id", "child_id")

        pathways$hsa_id <- map_pathways_to_ids(pathways, ptw_ids)
        pathways <- assign_hierarchy_levels(pathways, ptw_hr)
        pathway_data <- generate_pathway_data(pathways, ptw_titles)
        pathways_summary <- summarize_pathways(pathways, pathway_data)
        return(pathways_summary)
}



#' Load Pathway IDs from File
#'
#' @param pathways_hsa_ids_file Character string specifying the path to the file
#'   containing human pathway IDs.
#'
#' @return A `data.table` containing the loaded pathway IDs.
#'
#' @import data.table
#' @export
#'
load_pathway_ids <- function(pathways_hsa_ids_file) {
    if (!file.exists(pathways_hsa_ids_file)) {
        stop("Error: File 'reactome_pathways_hsa_id.txt' is missing.")
    }
    return(fread(pathways_hsa_ids_file))
}



#' Load Pathway Hierarchy
#'
#' Loads the pathway hierarchy from a specified file.
#'
#' @param pathways_hierachy_file Character. Path to the hierarchy file.
#'
#' @return A data.table containing the pathway hierarchy.
#' @export
load_pathway_hierarchy <- function(pathways_hierachy_file) {
  if (!file.exists(pathways_hierachy_file)) {
    stop("Error: File 'pathways_hierarchy.txt' is missing.")
  }
  return(fread(pathways_hierachy_file, header = FALSE))
}

#' Load Pathway Titles
#'
#' Loads the list of pathway titles from a specified file.
#'
#' @param list_of_pathways_file Character. Path to the pathways list file.
#'
#' @return A data.table containing the list of pathway titles.
#' @export
load_pathway_titles <- function(list_of_pathways_file) {
  if (!file.exists(list_of_pathways_file)) {
    stop("Error: File 'list_of_pathways.txt' is missing.")
  }
  return(fread(list_of_pathways_file, header = FALSE))
}


#' Map Pathways to HSA IDs
#'
#' Maps pathway names to their corresponding HSA IDs.
#'
#' @param pathways A data.frame or data.table with a 'pathway' column.
#' @param ptw_ids A data.frame or data.table with 'ptw_name_reactome' and
#'   'hsa_id' columns.
#'
#' @return A character vector of HSA IDs.
#' @export
map_pathways_to_ids <- function(pathways, ptw_ids) {
  hsa_ids <- ptw_ids$hsa_id[
    match(pathways$pathway, ptw_ids$ptw_name_reactome)
  ]
  if (any(is.na(hsa_ids))) {
    stop("Error: Some pathways could not be mapped to HSA IDs.")
  }
  return(hsa_ids)
}

#' Assign Hierarchy Levels
assign_hierarchy_levels <- function(pathways, ptw_hr) {
        current_level <- pathways$hsa_id
        for (i in 1:11) {
                next_level <- paste0("parent_id_level", i)
                pathways[[next_level]] <- ptw_hr$parent_id[
                        match(current_level, ptw_hr$child_id)
                ]
                current_level <- pathways[[next_level]]
        }
        return(pathways)
}

#' Generate Pathway Data
#'
#' Generates pathway metadata by linking pathways to their parent titles.
#'
#' @param pathways A data.table with columns 'pathway', 'hsa_id', and
#'   additional parent hierarchy columns.
#' @param ptw_titles A data.table with columns 'V1' (HSA ID) and 'V2'
#'   (pathway title).
#'
#' @return A data.table with columns: id, hsa_id, parent_hsa_id, and
#'   parent_ptw_title.
#' @export
#' 
generate_pathway_data <- function(pathways, ptw_titles) {
  ptws_dat <- NULL
  for (j in 1:nrow(pathways)) {
    ptw <- pathways[j, 2:ncol(pathways)]
    ptw <- unlist(ptw)
    ptw <- ptw[!is.na(ptw)]
    ptw <- ptw[length(ptw)]
    ptw_name <- ptw_titles$V2[match(ptw, ptw_titles$V1)]
    ptw_dat <- data.table(
      id = pathways$pathway[j],
      hsa_id = pathways$hsa_id[j],
      parent_hsa_id = ptw,
      parent_ptw_title = ptw_name
    )
    ptws_dat <- rbind(ptw_dat, ptws_dat)
  }
  return(ptws_dat)
}

#' Summarize Pathways
#'
#' Adds a category column to pathways using parent pathway titles.
#'
#' @param pathways A data.table with a 'pathway' column.
#' @param pathway_data A data.table with 'id' and 'parent_ptw_title' columns.
#'
#' @return The input pathways data.table with an added 'category' column.
#' @export
summarize_pathways <- function(pathways, pathway_data) {
  pathways$category <- pathway_data$parent_ptw_title[
    match(pathways$pathway, pathway_data$id)
  ]
  return(pathways)
}