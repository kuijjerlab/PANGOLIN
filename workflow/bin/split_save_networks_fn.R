library(gtools)
require(purrr)
require(tidyverse)
#' Read and combine LIONESS network files
#'
#' Reads individual .txt network files within a specified range and combines
#' them into a single data frame with samples as columns.
#'
#' @param networks_cancer Data frame with \code{file_path} and \code{tcga_id}
#'   columns containing network file paths and sample identifiers
#' @param start Integer. First network index to read
#' @param end Integer. Last network index to read
#'
#' @return Data frame with combined networks where rows are TF-gene edges
#'   and columns are samples (named by TCGA IDs)
#'
#' @export

read_networks_txt <- function(networks_cancer, start, end) {
    datalist <- 
        lapply(networks_cancer$file_path[start:end], function(x) fread(x))
    datalist <- lapply(datalist, function(x) melt(x))
    datalist <- map(datalist, ~ (.x %>% select(-variable)))
    net <- dplyr::bind_cols(datalist)
    colnames(net) <- networks_cancer$tcga_id[start:end]
    return(net)
}

#' Combine and save LIONESS networks for a cancer type
#'
#' Combines individual LIONESS network files for a specific cancer and saves
#' as RData. Uses memory-efficient splitting for large datasets (>=300 samples).
#'
#' @param tumor Character. Cancer type identifier (e.g., "BRCA", "LUAD")
#' @param info_net Data frame with \code{cancer}, \code{file_path}, and 
#'   \code{tcga_id} columns containing network metadata
#' @param outputdir Character. Directory path where RData files will be saved
#'
#' @return Invisible. Creates RData files containing combined networks
#'
#' @export

save_combined_networks <- function(tumor, info_net, outputdir) {
    cat("looking at", tumor, "\n")
    networks_cancer <- info_net %>% filter(cancer == tumor)
    if (nrow(networks_cancer)<300) {
        net <- read_networks_txt(networks_cancer, 1, nrow(networks_cancer))
        cat("Saving networks", "\n")
        save(net, file = file.path(outputdir, paste0("net_", tumor, ".RData")))
        rm(net); gc()
        } else {
            split_point_1 <- round(length(networks_cancer$file_path)/4)
            split_point_2 <- split_point_1*2
            split_point_3 <- split_point_1*3
            end_point <- length(networks_cancer$file_path)
            net <- read_networks_txt(networks_cancer, 1, split_point_1)
            cat("Saving the first part of networks", "\n")
            save(net,
                 file = file.path(outputdir, paste0("net_", tumor, "1.RData")))
            rm(net); gc()
            m <- split_point_1+1
            net <- read_networks_txt(networks_cancer, m, split_point_2)
            cat("Saving the second part of networks", "\n")
            save(net, 
                file = file.path(outputdir, paste0("net_", tumor, "2.RData")))
            rm(net); gc()
            d <- split_point_2+1
            cat("Saving the third part of networks")
            net <- read_networks_txt(networks_cancer, d, split_point_3)
            save(net, 
                file = file.path(outputdir, paste0("net_", tumor, "3.RData")))
            rm(net); gc()
            s <- split_point_3+1
            cat("Saving the forth part of networks", "\n")
            net <- read_networks_txt(networks_cancer, s, end_point)
            save(net, 
                file = file.path(outputdir, paste0("net_", tumor, "4.RData")))
        }
}