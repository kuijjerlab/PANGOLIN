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
#' This function runs the PORCUPINE methods on networks and pathway data for
#' a specific cancer type. It reads the network and edges files, filters the
#' pathways, performs PCA on the pathways, and runs permutations to assess
#' significance.
#'
#' @param cancer the type of cancer for which the analysis is to be performed.
#' @param network_dir the directory containing the network files.
#' @param edge_file A character string specifying the path to the file
#'                  containing edges information.
#' @param pathway_file A character string specifying the path to the file
#'                     containing pathway information in GMT format.
#' @param res_dir the directory where the results will be saved.
#' @param ncores_to_use the number of cores to use for parallel processing.
#' @return results are saved in a specified directory
#' @examples

runPORCUPINE_norm <- function(network_files,
                        cancer,
                        edge_file,
                        pathway_file,
                        res_dir,
                        ncores_to_use,
                        minSize = 5,
                        maxSize = 150) {
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
                file.path(res_dir, paste0("pathways_results_", cancer,
                                          ".txt")),
                col.names = T, row.names = F, sep = "\t", quote = FALSE)
    message("Running PORCUPINE permutations", "\n")
    pca_res_random <- pca_random(net, edges, pca_res_pathways, pathways,
                                 n_perm = 1000, ncores = ncores_to_use,
                                 scale_data = TRUE)

    write.table(pca_res_random,
                file.path(res_dir, paste0("pathways_results_random_",
                                          cancer, ".txt")),
                col.names = T, row.names = F, sep = "\t", quote = FALSE)
    res_porcupine <- porcupine(pca_res_pathways, pca_res_random)
    res_porcupine$p.adjust <- p.adjust(res_porcupine$pval, method = "fdr")
    write.table(res_porcupine,
                file.path(res_dir, paste0("pcp_results_", cancer, ".txt")),
                col.names = T, row.names = F, sep = "\t", quote = FALSE)
}

#' Read Sample File
#'
#' This function reads a sample file
#'
#' @param sample_file the path to the sample file to be read.
#' @return A data.table object containing the samples read from the file.
#' @export
#' 
read_sample_file <- function(sample_file) {
    samples <- fread(sample_file)
    return(samples)
}


#' Extract PD-1 Related Edges and Network Data
#'
#' This function extracts edges and network data related to PD-1 pathway
#'
#' @param cancer the type of cancer for which the analysis is to be performed.
#' @param network_dir the directory containing the network files.
#' @param edge_file the path to the file containing edges information. The 
#'                  file should have columns "reg" and "tar".
#' @param pathway_file the path to the file containing pathway information in
#'                     GMT format.
#' @param sample_file the path to a file containing sample information.
#' @param res_dir the directory where the extracted data will be saved.
#' @return saves the extracted PD-1 related network and edges data to 
#'         specified files in the result directory.
#' @export

extract_pd1_edges_norm <- function(network_files,
                        cancer,
                        edge_file,
                        pathway_file,
                        sample_file,
                        res_dir) {
    samples <- read_sample_file(sample_file)
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
    save(pd1_net, file = file.path(res_dir, paste0("pd1_net_norm_", cancer,
                                                   ".RData")))
    save(links, file = file.path(res_dir, paste0("pd1_edges_norm_", cancer,
                                                 ".RData")))
    }



#' Calculate Indegree for Each Cancer
#'
#' This function calculates the indegree for each cancer by reading network 
#' files and edge data.
#'
#' @param cancer the type of cancer for which the analysis is to be performed.
#' @param network_dir The directory where the network files are located.
#' @param edge_file The path to the txt file containing edges, with columns
#'                  "reg" and "tar".
#' @return Data frame. A table with the calculated indegrees for each
#'         target.
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



#' Extract mapping of individuals on PC1 for a pathway
#' @param cancer cancer 
#' @param network_dir Directory where the networks files are
#' @param edge_file Txt file with edges, columns "reg", "tar"
#' @param pathway_file Gmt file with pathways
#' @param sample_file File with samples
#' @param res_dir directory to store pcp results
#' @return Saving the ouput of pcp
#' @export

extract_individul_mapping_norm <- function(network_files,
                        cancer,
                        edge_file,
                        pathway_file,
                        sample_dir,
                        res_dir) {
    sample_files <- list.files(sample_dir)
    sample_file <- sample_files[grep(cancer, sample_files,
                                     ignore.case = TRUE)]
    load(file.path(sample_dir, sample_file), samples <- new.env())
    samples <- samples[["samples"]]
    message("Reading in network files", "\n")
    net <- read_networks_norm(network_files)
    colnames(net) <- samples$id
    message("Reading in edges file", "\n")
    edges <- fread(edge_file)
    pathways <- load_gmt(pathway_file)
    pathways <- filter_pathways(pathways, edges)
    pathways_to_use <- filter_pathways_size(pathways, minSize = 5,
                                             maxSize = 150)
    ind_scores <- get_pathway_ind_scores(pathways,
                        net,
                        edges)
    save(ind_scores, file = file.path(res_dir, paste0("individual_scores_",
                                                       cancer, ".RData")))
}



#' Extract mapping of individuals on PC1 for a pathway
#' @param cancer cancer 
#' @param network_dir Directory where the networks files are
#' @param edge_file Txt file with edges, columns "reg", "tar"
#' @param pathway_file Gmt file with pathways
#' @param sample_file File with samples
#' @param res_dir directory to store pcp results
#' @return Saving the ouput of pcp
#' @export

extract_feature_scores_norm <- function(network_files,
                        cancer,
                        edge_file,
                        pathway_file,
                        sample_dir,
                        res_dir) {
    sample_files <- list.files(sample_dir)
    sample_file <- sample_files[grep(cancer, sample_files,
                                     ignore.case = TRUE)]
    load(file.path(sample_dir, sample_file), samples <- new.env())
    samples <- samples[["samples"]]
    message("Reading in network files", "\n")
    net <- read_networks_norm(network_files)
    colnames(net) <- samples$id
    message("Reading in edges file", "\n")
    edges <- fread(edge_file)
    pathways <- load_gmt(pathway_file)
    pathways <- filter_pathways(pathways, edges)
    pathways_to_use <- filter_pathways_size(pathways, minSize = 5,
                                             maxSize = 150)
    feature_scores <- get_pathway_features_scores(pathways,
                        net,
                        edges)
    save(feature_scores, file = file.path(res_dir, paste0("feature_scores_",
                                                           cancer, ".RData")))
}