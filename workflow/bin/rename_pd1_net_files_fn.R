rename_pd1_net_file <- function(indegree_file, pd1_net_file){
    load(indegree_file, indegree <- new.env())
    indegree <- indegree[["ind"]]
    load(pd1_net_file, net <- new.env())
    pd1_net <- net[["pd1_net"]]
    colnames(pd1_net) <- colnames(indegree)[-1]
    return(pd1_net)
}


