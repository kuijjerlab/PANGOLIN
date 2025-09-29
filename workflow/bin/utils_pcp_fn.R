#' Read .RData network files and Combine them
#'
#' This function reads individual .RData network files from specified file 
#' paths and combines them together if there is more than one .RData file. 
#'
#' @param network_files vector of file paths to the network .RData files
#' @return A combined network object containing all the networks from the
#'         specified .RData files.
#' @export
#'
#' 
read_networks <- function(network_files) {
    # Ensure we have file paths
    if (length(network_files) == 0) {
        stop("No network files provided")
    }
    # Check if files exist
    missing_files <- network_files[!file.exists(network_files)]
    if (length(missing_files) > 0) {
        stop(sprintf("The following files do not exist:\n%s", 
                    paste(missing_files, collapse = "\n")))
    }
    cat("Following files will be processed:", "\n")
    print(basename(network_files))
    networks_all <- NULL
    for (i in seq_along(network_files)) {
        cat("Loading file", basename(network_files[i]), "\n")
        load(network_files[i], data <- new.env())
        networks <- data[["net"]]
        networks_all <- cbind(networks_all, networks)
        rm(networks); gc()
    }
    return(networks_all)
}
#' Quantile Normalize Networks
#'
#' This function reads network files from specified file paths and applies
#' quantile normalization to the combined network data.
#'
#' @param network_files vector of file paths to the network .RData files
#' @return A matrix representing the quantile normalized network.
#' @export
#'
quantile_normalize_net <- function(network_files) {
    cat("Reading in networks", "\n")
    net <- read_networks(network_files)
    cat("Normalizing network", "\n")
    net_norm <- normalize.quantiles(as.matrix(net))
    colnames(net_norm) <- colnames(net)
    return(net_norm)
    rm(net)
    gc()
}

#' Read quantile normalized networks from specified file paths
#'
#' This function reads normalized network data files from specified file 
#' paths. It concatenates the normalized networks from all files into a single
#' matrix.
#'
#' @param network_files vector of file paths to the normalized network .RData
#'                      files
#' @return A matrix containing the concatenated normalized network data from
#'         all specified files.
#' @export
#' 
read_networks_norm <- function(network_files) {
    # Ensure we have file paths
    if (length(network_files) == 0) {
        stop("No network files provided")
    }
    
    # Check if files exist
    missing_files <- network_files[!file.exists(network_files)]
    if (length(missing_files) > 0) {
        stop(sprintf("The following files do not exist:\n%s", 
                    paste(missing_files, collapse = "\n")))
    }
    
    cat("Following normalized files will be processed:", "\n")
    print(basename(network_files))
    
    networks_all <- NULL
    for (i in seq_along(network_files)) {
        cat("Loading file", basename(network_files[i]), "\n")
        load(network_files[i], data <- new.env())
        networks <- data[["net_norm"]]
        networks_all <- cbind(networks_all, networks)
        rm(networks); gc()
    }
    return(networks_all)
}

#' Run the PORCUPINE method on quantile normalized networks
#'
#' This function runs the PORCUPINE methods on networks and pathway data for a
#' specific cancer type. It reads the network and edges files, filters the
#' pathways, performs PCA on the pathways, and runs permutations to assess
#' significance.
#'
#' @param cancer Cancer type for which the analysis is performed.
#' @param network_files Vector of paths to normalized network .RData files.
#' @param edge_file Path to file containing edges information.
#' @param pathway_file Path to file containing pathway information in GMT 
#'                     format.
#' @param ncores_to_use Number of cores to use for parallel processing.
#' @param pathways_results_file Path to save pathway PCA results.
#' @param pathways_results_random_file Path to save random permutation results.
#' @param porcupine_results_file Path to save final PORCUPINE results.
#' @param individual_scores_file Path to save individual pathway scores.
#' @param minSize Minimum pathway size (default: 5).
#' @param maxSize Maximum pathway size (default: 150).
#' @param n_permutations Number of permutations (default: 1000).
#' @return Results are saved to specified output files.
#' @examples
#' # Example usage (not run)
#' @export

runPORCUPINE_norm <- function(network_files,
                        cancer,
                        edge_file,
                        pathway_file,
                        ncores_to_use,
                        pathways_results_file,
                        pathways_results_random_file,
                        porcupine_results_file,
                        individual_scores_file,
                        minSize = 5,
                        maxSize = 150,
                        n_permutations = 1000) {
    message("Reading in network files", "\n")
    net <- read_networks_norm(network_files)
    message("Reading in edges file", "\n")
    edges <- fread(edge_file)
    pathways <- load_gmt(pathway_file)
    pathways <- filter_pathways(pathways, edges)
    pathways_to_use <- filter_pathways_size(pathways,
                    minSize = minSize, maxSize = maxSize)
    message("Running PORCUPINE on pathways", "\n")
    pca_res_pathways <- pca_pathway(pathways_to_use, 
            net, edges, ncores_to_use, scale_data = TRUE)
    write.table(pca_res_pathways,
                pathways_results_file,
                col.names = T, row.names = F, sep = "\t", quote = FALSE)

    message("Calculating individual scores", "\n")
    ind_scores <- get_pathway_ind_scores(pathways,
                        net,
                        edges)
    # Save individual scores
    save(ind_scores, file = individual_scores_file)

    message("Running PORCUPINE permutations", "\n")
    pca_res_random <- pca_random(net, edges, pca_res_pathways, pathways,
                                 n_perm = n_permutations, ncores = ncores_to_use,
                                 scale_data = TRUE)

    write.table(pca_res_random,
                pathways_results_random_file,
                col.names = T, row.names = F, sep = "\t", quote = FALSE)
    message("Calculating PORCUPINE results", "\n")
    res_porcupine <- porcupine(pca_res_pathways, pca_res_random)
    res_porcupine$p.adjust <- p.adjust(res_porcupine$pval, method = "fdr")
    write.table(res_porcupine,
                porcupine_results_file,
                col.names = T, row.names = F, sep = "\t", quote = FALSE)
}

#' Read Sample File
#'
#' This function reads a sample file.
#'
#' @param sample_file Path to the sample file to be read.
#' @return A data.table object containing the samples read from the file.
#' @export
#' 
read_sample_file <- function(sample_file) {
    samples <- fread(sample_file)
    return(samples)
}


#' Extract PD-1 Related Edges and Network Data
#'
#' This function extracts edges and network data related to PD-1 pathway.
#'
#' @param cancer Cancer type for which the analysis is performed.
#' @param network_files Vector of paths to normalized network .RData files.
#' @param edge_file Path to file containing edges information. The file should
#'                  have columns "reg" and "tar".
#' @param pathway_file Path to file containing pathway information in GMT 
#'                     format.
#' @param pd1_edges_file Path to save extracted PD-1 edges data.
#' @param pd1_net_file Path to save extracted PD-1 network data.
#' @return Saves the extracted PD-1 related network and edges data to 
#'         specified files.
#' @export

extract_pd1_edges_norm <- function(network_files,
                        cancer,
                        edge_file,
                        pathway_file,
                        pd1_edges_file,
                        pd1_net_file) {
    message("Reading in network files", "\n")
    net <- read_networks_norm(network_files)
    message("Reading in edges file", "\n")
    edges <- fread(edge_file)
    pathways <- load_gmt(pathway_file)
    pd1_genes <- pathways[grep("PD_1", names(pathways))]
    pd1_genes <- pd1_genes[[1]]
    idx <- which(edges$tar %in% pd1_genes)
    pd1_net <- net[idx,]
    links <- edges[idx,]
    save(pd1_net, file = pd1_net_file)
    save(links, file = pd1_edges_file)
    }

#' Calculate Indegree for Each Cancer
#'
#' This function calculates the indegree for each cancer by reading network 
#' files and edge data.
#'
#' @param network_files Vector of paths to normalized network .RData files.
#' @param edge_file Path to txt file containing edges, with columns "reg" and
#'                  "tar".
#' @return Data frame with calculated indegrees for each target.
#' @export
#'
#' @import data.table
#' @import dplyr
#'
calculate_indegree_norm <- function(network_files, 
                            edge_file) {
    net <- read_networks_norm(network_files)
    edges <- fread(edge_file)
    net_edges <- cbind(edges, net)
    cat("Calculating indegree", "\n")
    indegree <-  net_edges %>%
                select(-c("reg")) %>%
                group_by(tar) %>%
                summarise_all(funs(sum))
    return(indegree)
    rm(net, edges, net_edges)
    gc()
}




