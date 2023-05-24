library(tidyverse)
library(igraph)
library(ggraph)
library(gridExtra)
library(graphlayouts)
set.seed(990212) # for reproducibility

# ------------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------------

# data_files <- list.files(path = "data/repeated_trials", pattern = "*.csv", full.names = TRUE) # day 1 and day 2
# data_files_known <- list.files(path = "data/human_expert_eval", pattern = "*.csv", full.names = TRUE) # 3 known patients
data_files <- list.files(path = "data/", pattern = "*.csv", full.names = TRUE)

# Swedish version of data_to_adjacency
data_to_adjacency <- function(patient_path){
  #' Creates an adjacency matrix with weights as elements (set to 0 for NAs and on
  #' diagonal) given the "patient_path" file path of the raw data, and a data frame
  #' with the node weights. The adjacency matrix and data frame with node weights are
  #' both returned in a list, in that order. Access the adjacency matrix by calling
  #' "data_to_adjacency(patient_path)[[1]]" with your patient file path, and access
  #' the node weights by calling "data_to_adjacency(patient_path)[[2]]" in the same way.
  #' This version has the Swedish symptom names + short descriptions in the adjacency
  #' matrices, which could be used for plotting networks.

  # making adjacency matrix from data
  patient <- read_csv2(patient_path, show_col_types = FALSE)
  colnames(patient)[1] <- "patient_id" # renaming for simplicity
  colnames(patient)[2] <- "relevance" # renaming for simplicity

  ## removing unnecessary columns and adding translated symptoms
  adjacency <- patient[!is.na(patient$relevance),] %>%
    select(-c(starts_with("x"), `Modify?`, Unkn, starts_with("..."))) %>%
    filter_at(1, all_vars(!is.na(.)))

  ## extracting node weights and creating data frame of node weights
  node_weights <- adjacency %>%
    select(patient_id, `Painful?`) %>%
    tail(-1) %>% # not a weight, just a string "Freq" not needed
    mutate(`Painful?` = as.numeric(`Painful?`)) %>%
    rename("node" = patient_id)

  ## removing symptoms not chosen by patient and columns no longer needed
  adjacency <- adjacency %>%
    filter(relevance != 0) %>% # removing symptoms not chosen by patient (row)
    select_if(function(col) !all(col[1] == 0)) %>% # removing symptoms not chosen by patient (col)
    select(-c(relevance, `Painful?`))

  ## selecting the symptoms chosen by patient
  node_weights <- node_weights%>%
    inner_join(adjacency, by = c("node" = "patient_id")) %>%
    select(node, `Painful?`) %>%
    rename("node_weights" = `Painful?`)

  ## removing column and rows no longer needed and converting to numeric
  adjacency <- adjacency[-1,] %>% # not needed
    select(- patient_id) %>%
    mutate_all(as.numeric) %>%
    replace(is.na(.), 0) %>% # there are no edges if weights NA, so we can replace NA with 0
    as.data.frame()

  rownames(adjacency) <- colnames(adjacency) # set row names

  # directed network where rows indicate starting node and columns end node (data is
  # given as the transpose of this)
  adjacency <- t(adjacency)

  list(as.data.frame(adjacency), node_weights) %>%
    return()
}

# English version
data_to_adjacency <- function(patient_path){
  #' Creates an adjacency matrix with weights as elements (set to 0 for NAs and on
  #' diagonal) given the "patient_path" file path of the raw data, and a data frame
  #' with the node weights. The adjacency matrix and data frame with node weights are
  #' both returned in a list, in that order. Access the adjacency matrix by calling
  #' "data_to_adjacency(patient_path)[[1]]" with your patient file path, and access
  #' the node weights by calling "data_to_adjacency(patient_path)[[2]]" in the same way.
  #' This version has the English symptom names + short descriptions in the adjacency
  #' matrices, which could be used for plotting networks.

  # making adjacency matrix from data
  patient <- read_csv2(patient_path, show_col_types = FALSE)
  colnames(patient)[1] <- "patient_id" # renaming for simplicity
  colnames(patient)[2] <- "relevance" # renaming for simplicity

  # vector of possible symptoms in English
  symptoms <- c("Insomnia", "Avoidance", "Feelings", "Overthinking", "Inactive", "Selfharm",
                "Substances", "Eating", "Suicidal", "Selfhate", "Unfocused", "Somatic",
                "School pressure", "Family situation", "Peer problems", "Trauma", "Dissociates", "Open item") %>%
    as_tibble()

  ## removing unnecessary columns and adding translated symptoms
  adjacency <- patient[!is.na(patient$relevance),] %>%
    select(- c(starts_with("x"), `Modify?`, Unkn, starts_with("..."))) %>%
    filter_at(1, all_vars(!is.na(.))) %>%
    mutate(symptoms = t(data_frame("symptoms", t(symptoms$value))))

  ## extracting node weights and creating data frame of node weights
  node_weights <- adjacency$`Painful?` %>%
    tail(-1) %>% # not a weight, just a string "Freq" not needed
    as.numeric()

  node_weights_df <- data_frame(t(data_frame(t(symptoms$value))), node_weights) # open item possible symptom
  names(node_weights_df)[1] = "node" # renaming for simplicity

  ## setting row and column names to translated symptoms and removing
  ## additional irrelevant columns
  adjacency <- adjacency %>%
    select(-`Painful?`) %>%
    relocate(symptoms) # move to first column

  adjacency <- adjacency %>%
    filter(relevance != 0) %>% # removing symptoms not chosen by patient (row)
    select_if(function(col) !all(col[1] == 0)) %>% # removing symptoms not chosen by patient (col)
    select(-patient_id, -relevance)

  adjacency <- adjacency[-1,] # not needed
  rownames <- t(adjacency[,1]) # store row names
  adjacency <- adjacency[,-1] # no longer needed

  adjacency <- adjacency %>% data.frame() # needed to be able to set row names

  rownames(adjacency) <- rownames # set row names
  colnames(adjacency) <- rownames(adjacency) # same column names as row names

  # removing symptoms not chosen by patient (col)
  adjacency <- adjacency %>%
    mutate_all(as.numeric) %>%
    replace(is.na(.), 0) %>% # there are no edges if weights NA, so we can replace NA with 0
    as.data.frame()

  # also removing unused node weights to keep only ones corresponding to relevant symptoms
  rownames(adjacency) <- rownames
  node_weights_df <- node_weights_df[which(node_weights_df$node %in% row.names(adjacency)),]

  # directed network where rows indicate starting node and columns end node (data is
  # given as the transpose of this)
  adjacency <- t(adjacency)

  list(as.data.frame(adjacency), node_weights_df) %>%
    return()
}

# ------------------------------------------------------------------------------
# Pruning methods
# ------------------------------------------------------------------------------

## ---------------------------------------------- ##
## Method 1: Naive approach with edge betweenness ##
## ---------------------------------------------- ##

pruning_edge_betweenness <- function(adjacency_node_weights, gamma = 0.95) {
  #' Computes the resulting graph when the method is applied to the network
  #' defined by the adjacency matrix and node weights in the list "adjacency_node_weight"
  #' (defined in this way to be used with the function "data_to_adjacency" above).
  #' The method uses the edge betweenness centrality according to Method I in the thesis.
  #' Can tune the value 0 <= "gamma" <= 1 to determine the level of pruning (0 means no
  #' pruning, 1 means that we have the same number of edges as the number of nodes).
  #' Returns a list of the pruned graph object and the original node weights, in that order.

  # extracting adjacency matrix and node weights
  adjacency <- adjacency_node_weights[[1]]
  node_weights_df <- adjacency_node_weights[[2]]

  graph <- graph_from_adjacency_matrix(as.matrix(adjacency),
                                       mode = "directed",
                                       weighted = TRUE)

  g <- graph # to save original "graph" and update "g" in loop below
  n <- ceiling(gamma * (length(E(g)) - length(V(g)))) # number of edges to prune, # edges = # nodes if gamma = 1
  n_edges <- length(E(g))
  edge_weights_keep <- E(g)$weight # used to

  if (n == 0) {
    # if there are no edges to prune, do nothing
    return(list(g, node_weights_df$node_weights))
  }

  # stopping when # nodes = # edges (> because we remove an edge within the loop)
  while (length(E(g)) > (n_edges - n)) {
    # transforming edge weights by incorporating node weights
    weights_df <- as_ids(head_of(g, E(g))) %>% # end nodes of each edge
      data_frame() %>%
      left_join(node_weights_df, by = c("." = "node")) %>% # finding end node weights

      # the end nodes are in the order of the nodes in E(g), so E(g)$weight gives the
      # corresponding weight of the edge connected to each end node, whose weight is
      # defined by node_weights

      # normalizing weights and computing new edge weights (NOTE - we divide by sum of E(graph)$weight, not E(g)$weight)
      mutate(new_edge_weights = node_weights / sum(node_weights) * E(g)$weight / sum(E(graph)$weight))

    # creating new weights and replacing the old weights in the graph object g
    new_edge_weights <- weights_df$new_edge_weights

    # edge betweenness with new weights (use 1 - weights because high weight = "short" path and
    # the function "edge_betweenness" interprets the weights as distances, and weights must be non-neg.)
    edge_betweenness_new_weights <- edge_betweenness(g, weights = 1 - new_edge_weights)

    # creating table of all edges in new_g with corresponding indices
    edge_list <- E(g) %>%
      as_ids() %>%
      as_tibble() %>%
      mutate(ind = seq_along(as_ids(E(g))))

    # saving the indices of edges to keep by ordering edge betweenness scores
    best_edges <- order(edge_betweenness_new_weights, decreasing = TRUE) %>%
      head(- 1) # remove worst edge (could have same edge betweenness?)

    removed_edge_ind <- order(edge_betweenness_new_weights, decreasing = TRUE) %>%
      tail(1)

    # new pruned subgraph with fewer edges
    subgraph <- subgraph.edges(g,
                               eids = edge_list$ind[-removed_edge_ind],
                               delete.vertices = FALSE) # update g

    # update starting graph and corresponding weights
    g <- subgraph
    edge_weights_keep <- new_edge_weights[edge_list$ind[-removed_edge_ind]]
  }

  E(g)$weight <- edge_weights_keep # updating weights of g

  return(list(g, node_weights_df$node_weights))
}

## ----------------------------------------------------------------------- ##
## Method 2: Naive approach with page rank first and then edge betweenness ##
## ----------------------------------------------------------------------- ##

pruning_page_rank <- function(adjacency_node_weights, gamma = 0.95) {
  #' Computes the resulting graph when the method is applied to the network found in
  #' the file whose file path is "data_file". The method uses the page rank values to
  #' update the node weights, according to Method II in the thesis. Can tune the value
  #' 0 <= "gamma" <= 1 to determine the level of pruning (0 means no pruning, 1 means
  #' that we have the same number of edges as the number of nodes). Returns a list of
  #' the graph object and the new page rank node weights, in that order.

  # extracting adjacency matrix and node weights
  adjacency <- adjacency_node_weights[[1]]
  node_weights_df <- adjacency_node_weights[[2]]

  starting_node_weights <- node_weights_df$node_weights

  g <- graph_from_adjacency_matrix(as.matrix(adjacency),
                                   mode = "directed",
                                   weighted = TRUE)

  # updating PageRanks and removing worst edge until # edges = # nodes
  new_g <- g

  n <- ceiling(gamma * (length(E(g)) - length(V(g)))) # number of edges to prune, # edges = # nodes
  n_edges <- length(E(new_g))
  edge_weights_keep <- E(g)$weight

  if (n == 0) {
    # if there are no edges to prune, do nothing
    return(list(g, node_weights_df$node_weights))
  }

  # default damping factor gamma = 0.85 and default theta = 1.
  nodes_pr <- page_rank(new_g, personalized = starting_node_weights, algo = "prpack")$vector

  # adding PageRank scores to the dataframe with node weights to get node names
  node_info_df <- nodes_pr %>%
    as.numeric() %>%
    data_frame(node_weights_df)
  names(node_info_df)[1] <- "page_rank" # renaming for simplicity

  while (length(E(new_g)) > (n_edges - n)) {
    edge_weight <- E(new_g)$weight / sum(E(g)$weight) # normalized (NOTE - we divide by sum of E(g)$weight, not E(new_g)$weight)

    # computing new edge weights based on PageRank scores instead of initial node weights
    weights_df <- as_ids(head_of(new_g, E(new_g))) %>% # end nodes
      data_frame() %>%
      left_join(node_info_df, by = c("." = "node")) %>%
      mutate(new_edge_weights = page_rank * edge_weight) # sum of PageRank is 1, so no need to normalize that

    # creating new weights and replacing the old weights in the graph object new_g
    new_edge_weights <- weights_df$new_edge_weights

    # edge betweenness with 1 - new edge weights because weights are interpreted as distances
    edge_betweenness_new_weights <- edge_betweenness(new_g, weights = 1 - new_edge_weights)
    # computing the indices of original edge_betweenness_new_weights in descending order
    # of edge betweenness and keeping the best ones

    edge_list <- E(new_g) %>%
      as_ids() %>%
      as_tibble() %>%
      mutate(ind = seq_along(as_ids(E(new_g)))) # all edges in new_g with corresp. inds.

    # saving the indices of edges to keep by ordering edge betweenness scores
    best_edges_no_start_ind <- order(edge_betweenness_new_weights,
                                     decreasing = TRUE) %>%
      head(- 1) # remove worst edge (could have same edge betweenness?)

    removed_edge_ind <- order(edge_betweenness_new_weights,
                              decreasing = TRUE) %>%
      tail(1)

    # new subgraph with fewer edges
    subgraph <- subgraph.edges(new_g,
                               eids = edge_list$ind[-removed_edge_ind],
                               delete.vertices = FALSE) # update g
    new_g <- subgraph # update starting graph
    edge_weights_keep <- new_edge_weights[edge_list$ind[-removed_edge_ind]]
  }

  E(new_g)$weight <- edge_weights_keep # updating weights of g

  return(list(new_g, node_info_df$page_rank))
}

## -------------------------------------------------------------------------- ##
## Method 3: Updating approach with page rank first and then edge betweenness ##
## -------------------------------------------------------------------------- ##

pruning_updating_page_rank <- function(adjacency_node_weights, gamma = 0.95){
  #' Computes the resulting graph when the method is applied to the network
  #' found in the file whose file path is "data_file". The method uses the page rank
  #' centrality to give new node weights according to Method III in the thesis.
  #' Can tune the value 0 <= "gamma" <= 1 to determine the level of pruning
  #' (0 means no pruning, 1 means that we have the same number of edges as the
  #' number of nodes). Returns a list of the graph object and the new page rank
  #' node weights, in that order.

  # extracting adjacency matrix and node weights
  adjacency <- adjacency_node_weights[[1]]
  node_weights_df <- adjacency_node_weights[[2]]

  starting_node_weights <- node_weights_df$node_weights

  g <- graph_from_adjacency_matrix(as.matrix(adjacency),
                                   mode = "directed",
                                   weighted = TRUE)

  # updating PageRanks and removing worst edge until # edges = # nodes
  new_g <- g

  n <- ceiling(gamma * (length(E(g)) - length(V(g)))) # number of edges to prune, # edges = # nodes
  n_edges <- length(E(g))
  edge_weights_keep <- E(g)$weight

  if (n == 0) {
    # if there are no edges to prune, do nothing
    return(list(g, node_weights_df$node_weights))
  }

  while (length(E(new_g)) > (n_edges - n)) {
    edge_weight <- E(new_g)$weight / sum(E(g)$weight)  # normalized (NOTE - we divide by sum of E(g)$weight, not E(new_g)$weight)

    # update PageRank in each step
    ## default damping factor gamma = 0.85 and default theta = 1.
    nodes_pr <- page_rank(new_g, personalized = starting_node_weights, algo = "prpack")$vector

    # adding PageRank scores to df with node weights to get node names
    node_info_df <- nodes_pr %>%
      as.numeric() %>%
      data_frame(node_weights_df)
    names(node_info_df)[1] <- "page_rank" # renaming for simplicity

    # computing new edge weights based on PageRank scores instead of initial node weights
    weights_df <- as_ids(head_of(new_g, E(new_g))) %>% # end nodes
      data_frame() %>%
      left_join(node_info_df, by = c("." = "node")) %>%
      mutate(new_edge_weights = page_rank * edge_weight) # sum of PageRank is 1, so no need to normalize that

    # creating new weights and replacing the old weights in the graph object new_g
    new_edge_weights <- weights_df$new_edge_weights

    # edge betweenness with 1 - new edge weights because weights are interpreted as distances
    edge_betweenness_new_weights <- edge_betweenness(new_g, weights = 1 - new_edge_weights)
    # computing the indices of original edge_betweenness_new_weights in descending order
    # of edge betweenness and keeping the best ones

    # finding the remaining edges in new_g that do not create start nodes
    edge_list <- E(new_g) %>%
      as_ids() %>%
      as_tibble() %>%
      mutate(ind = seq_along(as_ids(E(new_g)))) # all edges in new_g with corresp. inds.

    # saving the indices of edges to keep by ordering edge betweenness scores
    # for those edges that are known to not have created start nodes
    best_edges_no_start_ind <- order(edge_betweenness_new_weights
                                     decreasing = TRUE) %>%
      head(- 1) # remove worst edge (could have same edge betweenness?)

    # saving the index of the removed edge that is not known to create a start
    # node upon removal (the index of remaining_edges_ind)
    removed_edge_ind <- order(edge_betweenness_new_weights
                              decreasing = TRUE) %>%
      tail(1)

    # finding which of the original edges in new_g was removed
    removed_edge_ind <- remaining_edges_ind$ind[removed_edge_ind]

    removed_edge <- edge_list$value[removed_edge_ind]

    # new pruned subgraph with fewer edges
    subgraph <- subgraph.edges(new_g,
                               eids = edge_list$ind[-removed_edge_ind],
                               delete.vertices = FALSE) # update g
    new_g <- subgraph
    starting_node_weights <- node_info_df$page_rank # use previous PageRank values as new starting weights
    edge_weights_keep <- new_edge_weights[edge_list$ind[-removed_edge_ind]]
  }

  E(new_g)$weight <- edge_weights_keep # updating weights of g

  return(list(new_g, node_info_df$page_rank))
}

## ------------------------------------------------- ##
## Method 4: Brute force with connectivity of graph  ##
## ------------------------------------------------- ##

connectivity_nodes <- function(g, edge_weight) {
  #' This function computes the connectivity between each node in the input graph
  #' "g" defined as the sum of the inverse "length" of the shortest path between each pair
  #' of nodes, computed with respect to the edge weights. If the weights are 0 and
  #' produce infinite distances, as is the case on the diagonal, the inverse weights
  #' are set to 0 to be ignored in the sum computed in the function "connectivity_graph".

  # 1 - edge_weight since these are scaled between 0 and 1, and we want to minimize shortest paths, but
  # the edge_weights are interpreted as distances
  # 1 / distance since we want to minimize the distance, but maximize the connectivity. By removing edges,
  # the connectivity would decrease (or stay const.), while the distance would increase (or stay const.)
  mat <- 1 / distances(g, weights = 1 - edge_weight, mode = "out") # edges going out from nodes
  mat[!is.finite(mat)] <- 0 # infinite values are set to 0 to not impact the sum in "connectivity_graph"

  return(mat)
}

connectivity_graph <- function(g, edge_weight) {
  #' This function computes the connectivity of the graph "g" with given edge weights,
  #' using the function "connectivity_nodes" which computes the connectivity between
  #' each pair of nodes in the graph as a function of the shortest paths between the nodes.
  #' The graph is assumed to be directed and weighted, which is why the factor 2 is missing
  #' when comparing to the original brute force approach.

  (1 / (length(V(g)) * (length(V(g)) - 1))) * sum(connectivity_nodes(g, edge_weight)) %>%
    return()
}

pruning_connectivity_kept <- function(adjacency_node_weights, gamma = 0.95) {
  #' This function prunes the networks in a greedy fashion with the aim to keep
  #' the connectivity, according to Method IV in the thesis. It computes the
  #' resulting graph(s) of the file(s) defined in the list of files "data_file"
  #' and we can tune the value 0 <= "gamma" <= 1 to determine the level of pruning
  #' (0 means no pruning, 1 means that we have the same number of edges as the
  #' number of nodes). Returns a list of the graph object and the new page
  #' rank node weights, in that order.

  # extracting adjacency matrix and node weights
  adjacency <- adjacency_node_weights[[1]]
  node_weights_df <- adjacency_node_weights[[2]]

  graph <- graph_from_adjacency_matrix(as.matrix(adjacency),
                                       mode = "directed",
                                       weighted = TRUE)

  # 1 - edge weights because the function for connectivity uses 1 - edge weights,
  # so this cancels it out
  new_edge_weights <- 1 - E(graph)$weight

  if (unweighted == FALSE) {
    # transforming edge weights by incorporating node weights (could change relationship)
    weights_df <- as_ids(head_of(graph, E(graph))) %>% # end nodes of each edge
      data_frame() %>%
      # left_join(node_weights_df, by = c("." = "node")) %>% # finding end node weights
      left_join(node_weights_df, by = c("." = "node")) %>% # finding end node weights
      # the end nodes are in the order of the nodes in E(g), so E(g)$weight gives
      # corresponding weight of edge connected to each end node
      mutate(new_edge_weights = node_weights / sum(node_weights) * E(graph)$weight / sum(E(graph)$weight)) # normalize (NOTE - we divide by sum of E(graph)$weight, not E(g)$weight)
    new_edge_weights <- weights_df$new_edge_weights
  }

  g <- graph
  n <- ceiling(gamma * (length(E(g)) - length(V(g)))) # number of edges to prune, # edges = # nodes
  prod_rk <- 1
  n_pruned <- 0 # number of edges pruned

  if (n == 0) {
    # if there are no edges to prune, do nothing
    return(list(g, node_weights_df$node_weights))
  }

  # iteratively prune edge with highest rk value
  # for (r in 1:n) {
  while (n_pruned < n) {
    rk_largest <- - Inf
    edge_largest <- NULL

    for (i in 1:length(E(g))) {
      # removing an edge
      new_g <- delete_edges(g, as_ids(E(g)[i]))
      new_g_edge_weights <- new_edge_weights[-i]

      rk <- connectivity_graph(new_g, new_g_edge_weights) / connectivity_graph(g, new_edge_weights)

      # rk should be between 0 and 1
      if (rk > 1) {
        stop("rk is greater than 1")
      }

      if (rk > rk_largest) {
        # if removal of current edge preserves connectivity more than
        # previously tried edges
        rk_largest <- rk
        edge_largest <- as_ids(E(g)[i])
      }
    }

    # remove the edge with highest rk value
    j <- which(as_ids(E(g)) == edge_largest) # find index of which edge weight to remove
    new_g <- delete_edges(g, edge_largest)

    new_edge_weights_delete <- new_edge_weights[-j]
    prod_rk_delete <- prod_rk * rk_largest # connectivity kept pruned vs og graph

    n_pruned <- n_pruned + 1
    g <- new_g # updating g
    new_edge_weights <- new_edge_weights_delete # updating new_edge_weights
    prod_rk <- prod_rk_delete # updating prod_rk
  }
  # replacing old edge weights of g with the new updated edge weights
  if (unweighted == TRUE) {
    E(g)$weight <- 1 - new_edge_weights # changing back with 1 - (1 - E(g)$weight)
  }

  else {
    E(g)$weight <- new_edge_weights
  }

  return(list(g, node_weights_df$node_weights))
}

wrap_strings <- function(string_vec, width){
  #' Takes a vector of strings and the width of each line as arguments and wraps
  #' the string. Used in plots of networks.
  map_chr(string_vec, function(x) {
    paste(strwrap(x, width = width), collapse = "\n")
    }
  )
}

plot_weighted_network <- function(g, title, node_weights, edge_weights, edge_label = NULL, text_size = 2.5){
  #' Function for plotting the graph object "g" with the string "title" as the title
  #' of the plot, the nodes sized according to the numerical vector "node_weights",
  #' the edge widths according to the numerical vector "edge_weights", labels of
  #' edges "edge_label", and text size set by the numerical argument "text_size".
  #' Returns the resulting plot.

  size <- node_weights
  V(g)$weight <- node_weights # used for length of edges with arrows

  ## apply the function to wrap the node labels
  V(g)$name <- wrap_strings(V(g)$name, 12)

  # producing plot
  plot <- ggraph(g, layout = "stress") +
    ggtitle(title) +
    geom_edge_arc(color = "gray70",
                  strength = - 0.1,
                  arrow = arrow(angle = 15,
                                length = unit(0.15, "inches"),
                                ends = "last",
                                type = "closed"),
                  aes(edge_width = edge_weights,
                      start_cap = circle(node1.weight * 1.11, "pt"), # 1.2 to have the arrows start at node border
                      end_cap = circle(node2.weight * 1.2, "pt"), # 10 to get a space between end node and arrowhead
                      label = edge_label # possible edge label, default is NULL
                  ),
                  show.legend = c(edge_width = FALSE)) + # edge ending outside of node
    geom_node_point(color = "indianred1",
                    alpha = 0.8, # lower opacity to see edges behind nodes connected to different nodes
                    show.legend = c(size = FALSE, color = FALSE),
                    size = size) +
    geom_node_text(alpha = 1,
                   aes(label = as_ids(V(g))),
                   # aes(label = paste(as_ids(V(g)),"\n", round(node_weights * 100 / 30, 2))), # used if want label
                   repel = FALSE,
                   fontface = "bold",
                   colour = "black",
                   size = text_size
    ) +
    scale_edge_width(range = c(0, 1.5)) + # control size
    theme_graph(base_family = "serif", # otherwise error message about font
                # plot_margin = margin(70,70,70,70) # used for pdfs, not for other plots
    ) +
    coord_cartesian(clip = "off")
  return(plot)
}

plot_node_w_only_network <- function(g, title = "", node_weights = 100){
  #' Function for plotting the graph object "g" with the string "title" as the title
  #' of the plot, and the nodes sized according to the numerical vector "node_weights".
  #' Returns the resulting plot.

  size <- node_weights * 0.3
  V(g)$weight <- node_weights * 0.3 # used for length of edges with arrows

  ## apply the function to wrap the node labels
  V(g)$name <- wrap_strings(V(g)$name, 12)

  # producing plot
  plot <- ggraph(g, layout = "stress") +
    ggtitle(title) +
    geom_edge_arc(color = "gray70",
                  strength = - 0.1,
                  arrow = arrow(angle = 15,
                                length = unit(0.15, "inches"),
                                ends = "last",
                                type = "closed"),
                  aes(start_cap = circle(node1.weight * 1.11, "pt"), # 1.2 to have the arrows start at node border
                      # end_cap = circle(node2.weight * 1.3, "pt")), # 1.3 to get a space between end node and arrowhead for some plots
                      end_cap = circle(node2.weight + 10, "pt")), # 10 to get a space between end node and arrowhead for some other plots
                  show.legend = c(edge_width = FALSE)) + # edge ending outside of node
    geom_node_point(color = "indianred1",
                    alpha = 0.8, # lower opacity to see edges behind nodes connected to different nodes
                    show.legend = c(size = FALSE, color = FALSE),
                    size = size) +
    geom_node_text(alpha = 1,
                   aes(label = as_ids(V(g))),
                   repel = FALSE,
                   fontface = "bold",
                   colour = "black",
                   size = 2.5
    ) +
    theme_graph(base_family = "serif", # otherwise error message about font
                # plot_margin = margin(70,70,70,70) # used for pdfs, not for other plots
    ) +
    coord_cartesian(clip = "off")
  return(plot)
}

pruned_to_pdf <- function(data_file, file_name, gamma = 0.95) {
  #' Produces a pdf of the resulting networks for each method and the original
  #' network all on the same page with the original network at the top left. The
  #' pdf is saved as "file_name" and the data files of the original data are given
  #' by the list of strings "data_files". The value of "gamma" is set to 0.95.

  # Producing the plots for each method

  ## original plot
  plot_original <- map(data_file, function(x) {
    adjacency_node_weights <- data_to_adjacency(x)
    adjacency <- adjacency_node_weights[[1]]
    node_weights_df <- adjacency_node_weights[[2]]

    original_g <- graph_from_adjacency_matrix(as.matrix(adjacency),
                                              mode = "directed",
                                              weighted = TRUE)
    plot_weighted_network(original_g,
                          "Originalnätverk",
                          node_weights_df$node_weights / 100 * 30,
                          E(original_g)$weight / 100
    )
  })

  ## brute force network
  plot_pruning_connectivity_kept <- map(data_file, function(x) {
    adjacency_node_weights <- data_to_adjacency(x)
    node_weights_df <- adjacency_node_weights[[2]]

    graph_pruning_connectivity_kept <- pruning_connectivity_kept(data_to_adjacency(x), unweighted, gamma)

    plot_node_w_only_network(graph_pruning_connectivity_kept[[1]],
                 node_weights = node_weights_df$node_weights
    )
  })

  ## updating PageRank network
  plot_pruning_updating_page_rank <- map(data_file, function(x) {
    adjacency_node_weights <- data_to_adjacency(x)
    node_weights_df <- adjacency_node_weights[[2]]

    graph_pruning_updating_page_rank <- pruning_updating_page_rank(data_to_adjacency(x), unweighted, gamma)

    plot_node_w_only_network(graph_pruning_updating_page_rank[[1]],
                 node_weights = node_weights_df$node_weights
    )
  })

  # PageRank network
  plot_pr <- map(data_file, function(x){
    adjacency_node_weights <- data_to_adjacency(x)
    node_weights_df <- adjacency_node_weights[[2]]

    graph_pruning_page_rank <- pruning_page_rank(data_to_adjacency(x), unweighted, gamma)

    plot_node_w_only_network(graph_pruning_page_rank[[1]],
                 node_weights = node_weights_df$node_weights
    )
  })

  # edge betweenness network
  plot_edge_b <- map(data_file, function(x) {
    adjacency_node_weights <- data_to_adjacency(x)
    node_weights_df <- adjacency_node_weights[[2]]

    graph_pruning_edge_betweenness <- pruning_edge_betweenness(data_to_adjacency(x), unweighted, gamma)

    plot_node_w_only_network(graph_pruning_edge_betweenness[[1]],
                 node_weights = node_weights_df$node_weights
    ) # rescale
  })

  # arranging above plots in desired order
  arranged_plots <- c()
  for(i in 1:length(data_file)) {
    name <- paste("plot", i, sep = "")
    arranged_plots <- c(arranged_plots, assign(name, c(plot_edge_b[i],
                                                       plot_pr[i],
                                                       plot_pruning_updating_page_rank[i],
                                                       plot_pruning_connectivity_kept[i],
                                                       plot_original[i]
    )))

  }

  # setting desired layout, being 2 + 2 + 2 (3 rows, 2 columns)
  layout <- rbind(c(5,NA),
                  c(1,2),
                  c(3,4))

  # saving the plots
  ggsave(filename = file_name,
         plot = marrangeGrob(arranged_plots,
                             nrow = 3,
                             ncol = 2,
                             layout_matrix = layout,
                             top = quote(paste("Sida", g, "av", npages, "\n", "Patient-ID:",
                                               patientIDs[g]))),
         width = 15, height = 15)
}

# Finding which graphs were from day 1 and 2
day_1 <- data_files[c(TRUE, FALSE)]
day_2 <- data_files[c(FALSE, TRUE)]

# Change "data_files_known", "day_1", or "day_2" into the list of strings of data file names you want to use
patientIDs <- c("A-2", "B-2", "C-2") # use for 3 known patients
pruned_to_pdf(data_files_known, "IG_LA_LT.pdf") # use for 3 known patients

patientIDs <- str_extract(day_1, "(?<=data/).*(?=.csv)") # use for day 1
pruned_to_pdf(day_1, "day_1.pdf") # use for day 1

patientIDs <- str_extract(day_2, "(?<=data/).*(?=.csv)") # use for day 2
pruned_to_pdf(day_2, "day_2.pdf") # use for day 2

# The function is quite slow to run for large numbers of networks (such as for "day_1" or "day_2")

# ------------------------------------------------------------------------------
# Example networks
# ------------------------------------------------------------------------------

## -------------------------------------------------------- ##
## the example graph used in the introduction of the thesis ##
## -------------------------------------------------------- ##

example_adj <- rbind(c(0,100,0,0,0,0,0),
                     c(0,0,10,70,40,0,0),
                     c(0,20,0,0,0,0,0),
                     c(0,0,0,0,0,100,0),
                     c(0,0,0,0,0,0,0),
                     c(0,86,0,0,0,0,0),
                     c(0,0,0,0,0,0,0)) %>%
  as.matrix()

example_node_weights <- c(30,90,88,50,100,50,70)
example_graph <- graph_from_adjacency_matrix(example_adj, mode = "directed", weighted = TRUE)
V(example_graph)$name <- V(example_graph)

plot_weighted_network(example_graph,
                      "",
                      example_node_weights * 20 / 100,
                      E(example_graph)$weight * 20 / 100,
                      text_size = 3.2)

# --------------------------------------------------------
# the example graph used in the Dijkstra's alg
# --------------------------------------------------------

dijkstra_adj <- rbind(c(0,12,30,0,0),
                      c(0,0,10,86,0),
                      c(0,0,40,30,15),
                      c(0,0,0,0,60),
                      c(0,0,0,0,0)
                      ) %>%
  as.matrix()

dijkstra_node_weights <- c(90,90,90,90,90)
dijsktra_graph <- graph_from_adjacency_matrix(dijkstra_adj, mode = "directed", weighted = TRUE)
V(dijsktra_graph)$name <- V(dijsktra_graph)

plot_weighted_network(dijsktra_graph, "", dijkstra_node_weights * 20 / 100, 30 * 20 / 100,
                      edge_label = round(E(dijsktra_graph)$weight, 2), text_size = 4) # remember use labels in function

# --------------------------------------------------------
# example of hairball network
# --------------------------------------------------------

hairball <- data_to_adjacency(data_files[21])
hairball_graph <- graph_from_adjacency_matrix(hairball[[1]] %>% as.matrix(), mode = "directed", weighted = TRUE)

plot_weighted_network(hairball_graph, "",  hairball[[2]]$node_weights / 100 * 30,  E(hairball_graph)$weight / 100)

# ------------------------------------------------------------------------------
# example graphs for high and low node, edge betweenness, out-degree, page_rank?
# ------------------------------------------------------------------------------

plot_node_w_only_network_centralities <- function(g, title = "", node_weights = 20,
                                      color_nodes = "gray70", color_edges = "gray70",
                                      edge_centrality = FALSE, color_scale = NULL){
  #' Function for plotting the graph object "g" with the string "title" as the title
  #' of the plot, with the nodes sized according to the numerical vector "node_weights",
  #' the nodes colored according to "color_nodes", and the edges colored according to
  #' "color_edges". The logical variable "edge_centrality" determines whether the
  #' centrality to be plotted is an edge centrality or node centrality, and the variable
  #' "color_scale" determines the colors for the different values.

  # Load
  library(RColorBrewer)

  size <- node_weights
  V(g)$weight <- node_weights # used for length of edges with arrows
  my_palette <- colorRampPalette(brewer.pal(name = "Reds", n = 9)[3:9])

  ## apply the function to wrap the node labels
  V(g)$name <- wrap_strings(V(g)$name, 12)

  # producing plot
  plot <- ggraph(g, layout = "stress") +
    ggtitle(title) +
    {if(edge_centrality == FALSE)
      geom_edge_arc(color = "gray70",
                    strength = - 0.1,
                    edge_width = 1.5,
                    aes(start_cap = circle(node1.weight * 1.11, "pt"), # 1.11 to have the arrows start at node border
                        end_cap = circle(node2.weight * 1.11, "pt")), # 1.11 to get a space between end node and arrowhead
                    show.legend = c(edge_width = FALSE))} + # edge ending outside of node
    {if(edge_centrality == TRUE)
      geom_edge_arc(
        strength = - 0.1,
        arrow = arrow(angle = 15,
                      length = unit(0.1, "inches"),
                      ends = "last",
                      type = "closed"),
        edge_width = 1.5,
        aes(color = E(betweenness_graph)$max_between,
            start_cap = circle(node1.weight * 1.11, "pt"), # 1.11 to have the arrows start at node border
            end_cap = circle(node2.weight * 1.11, "pt")), # 1.11 to get a space between end node and arrowhead
        show.legend = c(edge_width = FALSE))} + # edge ending outside of node
    {if(edge_centrality == FALSE)
      geom_node_point(aes(color = V(g)$max_between),
                      alpha = 0.8, # lower opacity to see edges behind nodes connected to different nodes
                      show.legend = c(size = FALSE, color = FALSE),
                      size = size)} +
    {if(edge_centrality == TRUE)
      geom_node_point(color = "gray80",
                      alpha = 0.8, # lower opacity to see edges behind nodes connected to different nodes
                      show.legend = c(size = FALSE, color = FALSE),
                      size = size)} +
    theme_graph(base_family = "serif", # otherwise error message about font
    ) +
    {if(edge_centrality == FALSE)scale_color_manual(values = c("gray80", "indianred1"))} +
    {if(edge_centrality == TRUE)scale_edge_color_manual(values = color_scale)} +
    theme(legend.position = "none") +
    coord_cartesian(clip = "off")
  return(plot)
}

# edge-node-edge bridge
betweenness_adj <- rbind(c(0,1,0,0,0,0,0,0,0,0,0),
                         c(0,0,1,1,1,1,0,0,0,0,0),
                         c(0,0,0,0,0,0,0,0,0,0,0),
                         c(0,0,0,0,0,0,0,0,0,0,0),
                         c(0,0,0,0,0,0,0,0,0,0,0),
                         c(0,0,0,0,0,0,0,0,0,0,0),
                         c(0,0,0,0,0,0,0,1,0,0,0),
                         c(1,0,0,0,0,0,0,0,0,0,0),
                         c(0,0,0,0,0,0,0,1,0,0,0),
                         c(0,0,0,0,0,0,0,1,0,0,0),
                         c(0,0,0,0,0,0,0,1,0,0,0)) %>%
  as.matrix()

# node_edge_node bridge
betweenness_adj <- rbind(c(0,1,1,1,1,0,0,0,0,0),
                         c(0,0,0,0,0,0,0,0,0,0),
                         c(0,0,0,0,0,0,0,0,0,0),
                         c(0,0,0,0,0,0,0,0,0,0),
                         c(0,0,0,0,0,0,0,0,0,0),
                         c(0,0,0,0,0,0,1,0,0,0),
                         c(1,0,0,0,0,0,0,0,0,0),
                         c(0,0,0,0,0,0,1,0,0,0),
                         c(0,0,0,0,0,0,1,0,0,0),
                         c(0,0,0,0,0,0,1,0,0,0)) %>%
  as.matrix()

# loop as bridge
betweenness_adj <- rbind(c(0,1,1,1,1,1,0,0,0,0),
                         c(0,0,0,0,0,0,1,0,0,0),
                         c(0,0,0,0,0,0,0,0,0,0),
                         c(0,0,0,0,0,0,0,0,0,0),
                         c(0,0,0,0,0,0,0,0,0,0),
                         c(0,0,0,0,0,0,0,0,0,0),
                         c(1,0,0,0,0,0,0,1,1,1),
                         c(0,0,0,0,0,0,0,0,0,0),
                         c(0,0,0,0,0,0,0,0,0,0),
                         c(0,0,0,0,0,0,0,0,0,0)) %>%
  as.matrix() %>%
  t()

betweenness_graph <- graph_from_adjacency_matrix(betweenness_adj, mode = "directed", weighted = NULL)
V(betweenness_graph)$name <- V(betweenness_graph)
V(betweenness_graph)$betweeen <- betweenness(betweenness_graph)
V(betweenness_graph)$max_between <- ifelse(V(betweenness_graph)$betweeen == max(betweenness(betweenness_graph)),
                                           T, F)

E(betweenness_graph)$betweeen <- edge_betweenness(betweenness_graph)
E(betweenness_graph)$max_between <- ifelse(E(betweenness_graph)$betweeen %in% sort(edge_betweenness(betweenness_graph),
                                                                                   decreasing = T)[1],
                                           "max", ifelse(E(betweenness_graph)$betweeen %in% sort(edge_betweenness(betweenness_graph),
                                                                                                 decreasing = T)[2],
                                                         "2ndmax", ifelse(E(betweenness_graph)$betweeen %in% sort(edge_betweenness(betweenness_graph),
                                                                                                                  decreasing = T)[3],
                                                                          "3rdmax", F)))

plot_node_w_only_network_centralities(betweenness_graph,
                          color_nodes = as.factor(betweenness(betweenness_graph)))

plot_node_w_only_network_centralities(betweenness_graph,
                          color_edges = as.factor(edge_betweenness(betweenness_graph)),
                          edge_centrality = TRUE, color_scale = c("brown","gray70","lightpink", "indianred1") %>% rev())
plot_node_w_only_network_centralities(betweenness_graph,
                          color_edges = as.factor(edge_betweenness(betweenness_graph)),
                          edge_centrality = TRUE, color_scale = c("gray70","indianred1"))

# ------------------------------------------------------------------------------
# Repeated trials
# ------------------------------------------------------------------------------

## ------------------- ##
## Jaccard coefficient ##
## ------------------- ##

##### Jaccard index #####

jaccard <- function(A, B) {
  # Computes the Jaccard index for two sets A and B, defined as the cardinality
  # of the intersection of the two sets divided by the cardinality of their union.

  intersection <- length(intersect(A, B))
  union <- length(A) + length(B) - intersection
  return (intersection / union)
}

jaccard_tables_intersect <- function(adjacency_node_weights_1, adjacency_node_weights_2) {
  #' Computes a table with the Jaccard index computed for the graphs represented by
  #' "adjacency_node_weights_1" and "adjacency_node_weights_2" for all four methods
  #' and the original networks. We only consider the intersecting nodes when computing
  #' the Jaccard coefficients.

  intersecting_nodes <- intersect(adjacency_node_weights_1[[2]]$node, adjacency_node_weights_2[[2]]$node)

  # naive edge betweenness
  pruning_edge_betweenness_old <- pruning_edge_betweenness(adjacency_node_weights_1)[[1]] %>%
    subgraph(intersecting_nodes) %>%
    E() %>%
    as_ids()
  pruning_edge_betweenness_new <- pruning_edge_betweenness(adjacency_node_weights_2)[[1]] %>%
    subgraph(intersecting_nodes) %>%
    E() %>%
    as_ids()

  # edge betweenness with PageRank
  pruning_page_rank_old <- pruning_page_rank(adjacency_node_weights_1)[[1]] %>%
    subgraph(intersecting_nodes) %>%
    E() %>%
    as_ids()
  pruning_page_rank_new <- pruning_page_rank(adjacency_node_weights_2)[[1]] %>%
    subgraph(intersecting_nodes) %>%
    E() %>%
    as_ids()

  # Updated PageRank
  pruning_updating_page_rank_old <- pruning_updating_page_rank(adjacency_node_weights_1)[[1]] %>%
    subgraph(intersecting_nodes) %>%
    E() %>%
    as_ids()
  pruning_updating_page_rank_new <- pruning_updating_page_rank(adjacency_node_weights_2)[[1]] %>%
    subgraph(intersecting_nodes) %>%
    E() %>%
    as_ids()

  # Brute force
  pruning_connectivity_kept_old <- pruning_connectivity_kept(adjacency_node_weights_1)[[1]] %>%
    subgraph(intersecting_nodes) %>%
    E() %>%
    as_ids()
  pruning_connectivity_kept_new <- pruning_connectivity_kept(adjacency_node_weights_2)[[1]] %>%
    subgraph(intersecting_nodes) %>%
    E() %>%
    as_ids()

  # Jaccard
  jaccard <- list("original" = jaccard(graph_from_adjacency_matrix(adjacency_node_weights_1[[1]] %>% as.matrix(),
                                                                   mode = "directed",
                                                                   weighted = TRUE) %>%
                                         subgraph(intersecting_nodes) %>%
                                         E() %>%
                                         as_ids(),
                                       graph_from_adjacency_matrix(adjacency_node_weights_2[[1]] %>% as.matrix(),
                                                                   mode = "directed",
                                                                   weighted = TRUE) %>%
                                         subgraph(intersecting_nodes) %>%
                                         E() %>%
                                         as_ids()),
                  "edge_b" = jaccard(pruning_edge_betweenness_old, pruning_edge_betweenness_new),
                  "page_rank" = jaccard(pruning_page_rank_old, pruning_page_rank_new),
                  "updated_pr" = jaccard(pruning_updating_page_rank_old, pruning_updating_page_rank_new),
                  "connectivity_kept" = jaccard(pruning_connectivity_kept_old, pruning_connectivity_kept_new)
  )

  ## creating a table
  library(data.table)
  jaccard_to_table <- map(jaccard, as.data.table)
  jaccard_table <- rbindlist(jaccard_to_table, fill = TRUE, idcol = TRUE) %>%
    rename("method" = ".id", "jaccard_g" = "V1") %>%
    mutate(jaccard_g = jaccard_g %>% signif(3))

  return(jaccard_table %>% data_frame())
}


# Finding which graphs were from day 1 and 2
day_1 <- data_files[c(TRUE, FALSE)]
day_2 <- data_files[c(FALSE, TRUE)]

# Creating a list of the data frames of Jaccard indices for each patient
jaccard_tables_list <- list()
for (i in 1:length(day_1)) {
  assign(paste0("df", i), jaccard_tables_intersect(data_to_adjacency(day_1[i]), data_to_adjacency(day_2[i])))
}

jaccard_tables_list <- list(df1,df2,df3,df4,df5,df6,df7,df8,df9,df10,df11,
                            df12,df13,df14,df15,df16,df17,df18,df19,df20,df21)

# Creating data frame of all Jaccard indices
jaccard_table_all <- plyr::join_all(jaccard_tables_list, by = "method")

# Changing column names to be able to find row means
colnames <- names(jaccard_table_all)
for (i in 2:ncol(jaccard_table_all)) {
  colnames[i] <- paste0(colnames[i], "_", i - 1)
}
colnames(jaccard_table_all) <- colnames

# Tidy version of Jaccard table
jaccard_table_all_tidy <- jaccard_table_all %>%
  pivot_longer(2:ncol(jaccard_table_all), values_to = "jaccard", names_to = "person") %>%
  pivot_wider(names_from = "method", values_from = "jaccard") %>%
  pivot_longer(2:6, values_to = "jaccard", names_to = "method")

library(ggplot2)
# Needed to change names of facets
methods <- list("Edge betweenness", "PageRank", "Updated PageRank", "Connectivity kept")

labeller <- function(variable, value){
  return(methods[value])
}

library(ggpmisc)
# Facet scatter plot of relative Jaccards vs original Jaccards
jaccard_table_all %>%
  pivot_longer(2:ncol(jaccard_table_all), values_to = "jaccard", names_to = "person") %>%
  pivot_wider(names_from = "method", values_from = "jaccard") %>%
  mutate(pruned_edge_b = edge_b,
         pruned_page_rank = page_rank,
         pruned_updated_page_rank = updated_pr,
         pruned_connectivity_kept = connectivity_kept) %>%
  # filter(person != "jaccard_g_6") %>% # outlier
  select(2, 7:10) %>%
  pivot_longer(2:5, values_to = "jaccard", names_to = "method") %>%
  mutate(method = factor(method,
                         levels = c("pruned_edge_b", "pruned_page_rank", "pruned_updated_page_rank", "pruned_connectivity_kept"))) %>%
  ggplot(aes(original,
             jaccard)) +
  geom_point() +
  geom_smooth(se = FALSE, method = "lm") +
  stat_poly_eq(aes(label =  ..rr.label..)) +
  xlim(c(0, 1)) +
  xlab("Original Jaccard coefficient") +
  ylab("Pruned Jaccard coefficient") +
  scale_y_continuous(breaks = seq(-0.2, 1, 0.2)) +
  scale_x_continuous(breaks = seq(-0.2, 1, 0.2)) +
  facet_wrap(~method,
             labeller = labeller
  ) +
  ggtitle("Original vs. pruned Jaccard coefficients for repeated trials")

# Computing median without outlier
jaccard_table_all_tidy %>%
  filter(!(person == "jaccard_g_9" & method == "connectivity_kept")) %>%
  filter(method == "connectivity_kept") %>%
  summarize(median(jaccard))

# Percentage and number of network pairs with higher Jaccard after pruning
jaccard_table_all %>%
  pivot_longer(2:22, values_to = "jaccard", names_to = "person") %>%
  pivot_wider(names_from = "method", values_from = "jaccard") %>%
  mutate(pruned_edge_b = edge_b / original,
         pruned_page_rank = page_rank / original,
         pruned_updated_page_rank = updated_pr / original,
         pruned_connectivity_kept = connectivity_kept / original) %>%
  select(1:2, 7:10) %>%
  pivot_longer(3:6, values_to = "jaccard", names_to = "method") %>%
  mutate(increase = ifelse(jaccard >= 1, TRUE,
                           ifelse(jaccard < 1, FALSE, NA))) %>%
  group_by(method, increase) %>%
  count(increase) %>%
  filter(increase == TRUE)

# Boxplot
jaccard_table_all_pruning <- jaccard_table_all %>%
  pivot_longer(-method ,names_to = "person", values_to = "jaccard") %>%
  filter(method != "original")

jaccard_table_all_original <- jaccard_table_all %>%
  pivot_longer(-method ,names_to = "person", values_to = "jaccard") %>%
  filter(method == "original")

jaccard_table_all_boxplot <- jaccard_table_all_original %>%
  full_join(jaccard_table_all_pruning, by = "person") %>%
  rename("method" = "method.y", "jaccard_pruned" = "jaccard.y", "jaccard_original" = "jaccard.x") %>%
  select(-1) %>%
  pivot_longer(c(jaccard_original, jaccard_pruned), names_to = "type", values_to = "jaccard")

jaccard_table_all_boxplot %>%
  ggplot(aes(factor(type),
             jaccard)) +
  geom_boxplot() +
  facet_wrap(~ factor(method,
                      levels = c("edge_b", "page_rank", "updated_pr", "connectivity_kept"),
                      labels = c("Edge betweenness", "PageRank", "Updated PageRank", "Connectivity kept")),
             ncol = 4) +
  scale_x_discrete(limits = c("jaccard_pruned", "jaccard_original"),
                   labels = c("Pruned", "Original")) +
  scale_y_continuous(breaks = seq(-0.2, 1, 0.2)) +
  xlab("") +
  ylab("Jaccard coefficient") +
  ggtitle("Distributions of Jaccard coefficients for repeated trials")


## -------------------- ##
## Pearson and Spearman ##
## -------------------- ##

correlation <- function(g_1, g_2, method = "pearson") {
  #' Computes the Pearson or Spearman correlation coefficients for the adjacency
  #' matrix for the given graph and the associated node weights, depending on the
  #' argument "method" which can be set to "pearson" or "spearman". The arguments
  #' g_1 and g_2 should be lists containing the the adjacency matrix in the first
  #' position and the node weights in the second position.

  adjacency_old <- g_1[[1]]
  adjacency_new <- g_2[[1]]

  node_weights_old <- g_1[[2]]
  node_weights_new <- g_2[[2]]

  # correlation
  adjacency <- cor(adjacency_old %>% c(), adjacency_new %>% c(), method = method)
  painful <- cor(node_weights_old, node_weights_new, method = method)

  list(adjacency, painful) %>%
    return()
}

correlation_ranked <- function(g_1, g_2, method = "pearson") {
  #' Computes the Spearman correlation coefficients (ranked Pearson) for the adjacency
  #' matrix for the given graph and the associated node weights. The arguments
  #' g_1 and g_2 should be lists containing the the adjacency matrix in the first
  #' position and the node weights in the second position. (Can also set the argument
  #' "method" = "spearman" for this, but this way we ensure we know how Spearman is
  #' computed when we have ties.)

  adjacency_old <- g_1[[1]]
  adjacency_new <- g_2[[1]]

  node_weights_old <- g_1[[2]]
  node_weights_new <- g_2[[2]]

  # correlation (ranked with average for ties)
  adjacency <- cor(adjacency_old %>% c() %>% rank(), adjacency_new %>% c() %>% rank(), method = method)
  painful <- cor(node_weights_old %>% rank(), node_weights_new %>% rank(), method = method)

  list(adjacency, painful) %>%
    return()
}

get_intersections <- function(adjacency_node_weights_1, adjacency_node_weights_2) {
  #' Gets the submatrices of the adjacency matrices in "adjacency_node_weights_1"
  #' and "adjacency_node_weights_2" such that they only contain the intersecting
  #' symptoms chosen for both cases. Returns the corresponding node weights in a
  #' list of lists to correspond to the previously defined functions.

  original_old_adj <- adjacency_node_weights_1[[1]] %>%
    as.matrix()
  original_old_weights <- adjacency_node_weights_1[[2]]

  original_new_adj <- adjacency_node_weights_2[[1]] %>%
    as.matrix()
  original_new_weights <- adjacency_node_weights_2[[2]]

  intersect_symptoms <- intersect(rownames(original_old_adj), rownames(original_new_adj))

  ## subsetting adjacency matrices and node weights
  original_old_adj <- original_old_adj[rownames(original_old_adj) %in% intersect_symptoms,
                                       colnames(original_old_adj) %in% intersect_symptoms]

  original_new_adj <- original_new_adj[rownames(original_new_adj) %in% intersect_symptoms,
                                       colnames(original_new_adj) %in% intersect_symptoms]

  original_old_weights <- original_old_weights %>% filter(node %in% intersect_symptoms)

  original_new_weights <- original_new_weights %>% filter(node %in% intersect_symptoms)

  return(list(list(original_old_adj, original_old_weights),
              list(original_new_adj, original_new_weights)))
}

correlation_tables_intersect <- function(adjacency_node_weights_1, adjacency_node_weights_2, method, corel, ...) {
  #' Computes a correlation table between "adjacency_node_weights_1" and "adjacency_node_weights_2"
  #' for each of the results of the pruning methods as well as the original data. The table is
  #' returned at the end and one can choose the correlation coefficient to be computed by setting
  #' "method" to "pearson" or "spearman" (but it should be "pearson" in our cases). The function "corel"
  #' defines which correlation coefficient should be computed (correlation or correlation_ranked).

  # computing original adjacency and node weights to find original correlation
  intersections <- get_intersections(adjacency_node_weights_1, adjacency_node_weights_2)

  original_old <- intersections[[1]]
  original_new <- intersections[[2]]

  ## subsetting adjacency matrices and node weights
  original_old_adj <- original_old[[1]]
  original_new_adj <- original_new[[1]]

  original_old_weights <- original_old[[2]]$node_weights
  original_new_weights <- original_new[[2]]$node_weights

  ## finding correlations using correlation function "cor"
  original_correlations <- corel(list(original_old_adj, original_old_weights),
                                 list(original_new_adj, original_new_weights),
                                 method)

  # naive edge betweenness
  pruning_edge_betweenness_old <- pruning_edge_betweenness(adjacency_node_weights_1)
  pruning_edge_betweenness_new <- pruning_edge_betweenness(adjacency_node_weights_2)
  intersections <- get_intersections(list(pruning_edge_betweenness_old[[1]] %>%
                                            as_adjacency_matrix(attr = "weight", sparse = FALSE) %>%
                                            as.data.frame(),
                                          data_frame(node_weights = pruning_edge_betweenness_old[[2]],
                                                     node = pruning_edge_betweenness_old[[1]] %>% V() %>% as_ids)),
                                     list(pruning_edge_betweenness_new[[1]] %>%
                                            as_adjacency_matrix(attr = "weight", sparse = FALSE) %>%
                                            as.data.frame(),
                                          data_frame(node_weights = pruning_edge_betweenness_new[[2]],
                                                     node = pruning_edge_betweenness_new[[1]] %>% V() %>% as_ids))
  )

  pruning_edge_betweenness_old_full <- intersections[[1]]
  pruning_edge_betweenness_new_full <- intersections[[2]]

  ## subsetting adjacency matrices and node weights
  pruning_edge_betweenness_old_adj <- pruning_edge_betweenness_old_full[[1]]
  pruning_edge_betweenness_new_adj <- pruning_edge_betweenness_new_full[[1]]

  pruning_edge_betweenness_old_weights <- pruning_edge_betweenness_old_full[[2]]$node_weights
  pruning_edge_betweenness_new_weights <- pruning_edge_betweenness_new_full[[2]]$node_weights

  # edge betweenness with PageRank
  pruning_page_rank_old <- pruning_page_rank(adjacency_node_weights_1)
  pruning_page_rank_new <- pruning_page_rank(adjacency_node_weights_2)
  intersections <- get_intersections(list(pruning_page_rank_old[[1]] %>%
                                            as_adjacency_matrix(attr = "weight", sparse = FALSE) %>%
                                            as.data.frame(),
                                          data_frame(node_weights = pruning_page_rank_old[[2]],
                                                     node = pruning_page_rank_old[[1]] %>% V() %>% as_ids)),
                                     list(pruning_page_rank_new[[1]] %>%
                                            as_adjacency_matrix(attr = "weight", sparse = FALSE) %>%
                                            as.data.frame(),
                                          data_frame(node_weights = pruning_page_rank_new[[2]],
                                                     node = pruning_page_rank_new[[1]] %>% V() %>% as_ids))
  )

  pruning_page_rank_old_full <- intersections[[1]]
  pruning_page_rank_new_full <- intersections[[2]]

  ## subsetting adjacency matrices and node weights
  pruning_page_rank_old_adj <- pruning_page_rank_old_full[[1]]
  pruning_page_rank_new_adj <- pruning_page_rank_new_full[[1]]

  pruning_page_rank_old_weights <- pruning_page_rank_old_full[[2]]$node_weights
  pruning_page_rank_new_weights <- pruning_page_rank_new_full[[2]]$node_weights

  # Updated PageRank
  pruning_updating_page_rank_old <- pruning_updating_page_rank(adjacency_node_weights_1)
  pruning_updating_page_rank_new <- pruning_updating_page_rank(adjacency_node_weights_2)
  intersections <- get_intersections(list(pruning_updating_page_rank_old[[1]] %>%
                                            as_adjacency_matrix(attr = "weight", sparse = FALSE) %>%
                                            as.data.frame(),
                                          data_frame(node_weights = pruning_updating_page_rank_old[[2]],
                                                     node = pruning_updating_page_rank_old[[1]] %>% V() %>% as_ids)),
                                     list(pruning_updating_page_rank_new[[1]] %>%
                                            as_adjacency_matrix(attr = "weight", sparse = FALSE) %>%
                                            as.data.frame(),
                                          data_frame(node_weights = pruning_updating_page_rank_new[[2]],
                                                     node = pruning_updating_page_rank_new[[1]] %>% V() %>% as_ids))
  )

  pruning_updating_page_rank_old_full <- intersections[[1]]
  pruning_updating_page_rank_new_full <- intersections[[2]]

  ## subsetting adjacency matrices and node weights
  pruning_updating_page_rank_old_adj <- pruning_updating_page_rank_old_full[[1]]
  pruning_updating_page_rank_new_adj <- pruning_updating_page_rank_new_full[[1]]

  pruning_updating_page_rank_old_weights <- pruning_updating_page_rank_old_full[[2]]$node_weights
  pruning_updating_page_rank_new_weights <- pruning_updating_page_rank_new_full[[2]]$node_weights

  # Brute force
  pruning_connectivity_kept_old <- pruning_connectivity_kept(adjacency_node_weights_1)
  pruning_connectivity_kept_new <- pruning_connectivity_kept(adjacency_node_weights_2)
  intersections <- get_intersections(list(pruning_connectivity_kept_old[[1]] %>%
                                            as_adjacency_matrix(attr = "weight", sparse = FALSE) %>%
                                            as.data.frame(),
                                          data_frame(node_weights = pruning_connectivity_kept_old[[2]],
                                                     node =  pruning_connectivity_kept_old[[1]] %>% V() %>% as_ids)),
                                     list(pruning_connectivity_kept_new[[1]] %>%
                                            as_adjacency_matrix(attr = "weight", sparse = FALSE) %>%
                                            as.data.frame(),
                                          data_frame(node_weights = pruning_connectivity_kept_new[[2]],
                                                     node = pruning_connectivity_kept_new[[1]] %>% V() %>% as_ids))
  )

  pruning_connectivity_kept_old_full <- intersections[[1]]
  pruning_connectivity_kept_new_full <- intersections[[2]]

  ## subsetting adjacency matrices and node weights
  pruning_connectivity_kept_old_adj <- pruning_connectivity_kept_old_full[[1]]
  pruning_connectivity_kept_new_adj <- pruning_connectivity_kept_new_full[[1]]

  pruning_connectivity_kept_old_weights <- pruning_connectivity_kept_old_full[[2]]$node_weights
  pruning_connectivity_kept_new_weights <- pruning_connectivity_kept_new_full[[2]]$node_weights

  # correlation
  correlation <- list("original" = original_correlations,
                      "edge_b" = corel(list(pruning_edge_betweenness_old_adj, pruning_edge_betweenness_old_weights),
                                       list(pruning_edge_betweenness_new_adj, pruning_edge_betweenness_new_weights),
                                       method),
                      "page_rank" = corel(list(pruning_page_rank_old_adj, pruning_page_rank_old_weights),
                                          list(pruning_page_rank_new_adj, pruning_page_rank_new_weights),
                                          method),
                      "updated_pr" = corel(list(pruning_updating_page_rank_old_adj, pruning_updating_page_rank_old_weights),
                                           list(pruning_updating_page_rank_new_adj, pruning_updating_page_rank_new_weights),
                                           method),
                      "connectivity_kept" = corel(list(pruning_connectivity_kept_old_adj, pruning_connectivity_kept_old_weights),
                                            list(pruning_connectivity_kept_new_adj, pruning_connectivity_kept_new_weights),
                                            method)
  )
  ## creating a table
  library(data.table)
  correlation_to_table <- map(correlation, as.data.table)
  correlation_table <- rbindlist(correlation_to_table, fill = TRUE, idcol = T) %>%
    rename("method" = ".id", "cor_g" = "V1", "cor_node_weights" = "V2") %>%
    mutate(cor_g = cor_g %>% signif(3),
           cor_node_weights = cor_node_weights %>% signif(3),
           rel_cor_g = (cor_g / cor_g[1]) %>% signif(3),
           rel_cor_node_weights = (cor_node_weights / cor_node_weights[1]) %>% signif(3))

  return(correlation_table)
}

# Creating a list of the data frames of cor indices for each patient
cor_tables <- function(f) {
  cor_tables_list <- list()
  for (i in 1:length(day_1)) {
    assign(paste0("df", i),
           correlation_tables_intersect(data_to_adjacency(day_1[i]),
                                        data_to_adjacency(day_2[i]),
                                        method = "pearson",
                                        f))
  }

  cor_tables_list <- list(df1,df2,df3,df4,df5,df6,df7,df8,df9,df10,df11,
                          df12,df13,df14,df15,df16,df17,df18,df19,df20,df21)

  # Creating data frame of all cor indices
  cor_table_all <- plyr::join_all(cor_tables_list, by = "method")

  # Changing column names to be able to find row means
  colnames <- names(cor_table_all)
  for (i in 2:ncol(cor_table_all)) {
    colnames[i] <- paste0(colnames[i], "_", i - 1)
  }
  colnames(cor_table_all) <- colnames

  cor_table_all_tidy <- cor_table_all %>%
    select(method, starts_with("cor_g")) %>%
    pivot_longer(2:22,values_to = "cor", names_to = "person") %>%
    pivot_wider(names_from = "method", values_from = "cor") %>%
    pivot_longer(2:6, values_to = "pearson", names_to = "method")

  return(list(cor_table_all, cor_table_all_tidy))
}

# the Pearson correlation results when not using ranks (intersecting symptoms only)
cor_table_pearson <- cor_tables(correlation)
cor_table_all <- cor_table_pearson[[1]]
cor_table_all_tidy <- cor_table_pearson[[2]]
## this person has no intersecting edges for updating pagerank, so the correlation is NA and we
## thus replace it with 0 as the matrices are not correlated at all if there is nothing to check
cor_table_all$cor_g_13[is.na(cor_table_all$cor_g_13)] <- 0
cor_table_all_tidy$pearson[is.na(cor_table_all_tidy$pearson)] <- 0

# the Spearman correlation results when using ranks (intersecting symptoms only)
cor_table_spearman <- cor_tables(correlation_ranked)
cor_table_all <- cor_table_spearman[[1]]
cor_table_all_tidy <- cor_table_spearman[[2]]
## this person has no intersecting edges for updating pagerank, so the correlation is NA and we
## thus replace it with 0 as the matrices are not correlated at all if there is nothing to check
cor_table_all$cor_g_13[is.na(cor_table_all$cor_g_13)] <- 0
cor_table_all_tidy$pearson[is.na(cor_table_all_tidy$pearson)] <- 0

# Pearson correlation boxplot

# Boxplot
cor_table_all_pruning <- cor_table_all %>%
  select(method, starts_with("cor_g")) %>%
  pivot_longer(-method ,names_to = "person", values_to = "correlation") %>%
  filter(method != "original")

cor_table_all_original <- cor_table_all %>%
  select(method, starts_with("cor_g")) %>%
  pivot_longer(-method ,names_to = "person", values_to = "correlation") %>%
  filter(method == "original")

cor_table_all_boxplot <- cor_table_all_original %>%
  full_join(cor_table_all_pruning, by = "person") %>%
  rename("method" = "method.y", "correlation_pruned" = "correlation.y", "correlation_original" = "correlation.x") %>%
  select(-1) %>%
  pivot_longer(c(correlation_original, correlation_pruned), names_to = "type", values_to = "correlation")

cor_table_all_boxplot %>%
  ggplot(aes(factor(type),
             correlation)) +
  geom_boxplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 0.65) +
  facet_wrap(~ factor(method,
                      levels = c("edge_b", "page_rank", "updated_pr", "connectivity_kept"),
                      labels = c("Edge betweenness", "PageRank", "Updated PageRank", "Connectivity kept")),
             ncol = 4) +
  scale_x_discrete(limits = c("correlation_pruned", "correlation_original"),
                   labels = c("Pruned", "Original")) +
  scale_y_continuous(breaks = seq(-0.2, 1, 0.2)) +
  xlab("") +
  ylab("Spearman coefficient") +
  ggtitle("Distributions of Spearman coefficients for repeated trials")

check_linear <- function(pruning_method) {
  qq <- cor_table_all %>%
    pivot_longer(2:ncol(cor_table_all), values_to = "pearson", names_to = "person") %>%
    pivot_wider(names_from = "method", values_from = "pearson") %>%
    mutate(pruned_edge_b = edge_b,
           pruned_page_rank = page_rank,
           pruned_updated_page_rank = updated_pr,
           pruned_connectivity_kept = connectivity_kept) %>%
    select(2, 7:10) %>%
    pivot_longer(2:5, values_to = "pearson", names_to = "method") %>%
    filter(method == pruning_method) %>%
    select(original, pearson)

  lm(qq, formula = original ~ pearson) %>% plot()
}
check_linear("pruned_edge_b")
check_linear("pruned_page_rank")
check_linear("pruned_updated_page_rank")
check_linear("pruned_connectivity_kept")


# Pearson correlation plot of relative change in comparison to original networks
cor_table_all %>%
  select(method, starts_with("cor_g")) %>%
  pivot_longer(2:22, values_to = "pearson", names_to = "person") %>%
  pivot_wider(names_from = "method", values_from = "pearson") %>%
  mutate(pruned_edge_b = edge_b,
         pruned_page_rank = page_rank,
         pruned_updated_page_rank = updated_pr ,
         pruned_connectivity_kept = connectivity_kept ) %>%
  select(2, 7:10) %>%
  pivot_longer(2:5, values_to = "pearson", names_to = "method") %>%
  mutate(method = factor(method,
                         levels = c("pruned_edge_b", "pruned_page_rank", "pruned_updated_page_rank", "pruned_connectivity_kept"))) %>%
  ggplot(aes(original,
             pearson)) +
  geom_point() +
  geom_smooth(se = FALSE, method = "lm") +
  stat_poly_eq(aes(label =  ..rr.label..)) +
  xlim(c(NA, 1)) +
  xlab("Original Spearman coefficient") +
  ylab("Pruned Spearman coefficient") +
  scale_y_continuous(breaks = seq(-0.2, 1, 0.2)) +
  scale_x_continuous(breaks = seq(-0.2, 1, 0.2)) +
  facet_wrap(~method,
             labeller = labeller
  ) +
  ggtitle("Original vs. pruned Spearman coefficients for repeated trials
          ")

# Pearson correlation percentage of increased correlations
cor_table_all %>%
  select(method, starts_with("cor_g")) %>%
  pivot_longer(2:22, values_to = "pearson", names_to = "person") %>%
  pivot_wider(names_from = "method", values_from = "pearson") %>%
  drop_na() %>%
  mutate(pruned_edge_b = edge_b >= original,
         pruned_page_rank = page_rank >= original,
         pruned_updated_page_rank = updated_pr >= original,
         pruned_connectivity_kept = connectivity_kept >= original) %>%
  select(1:2, 7:10) %>%
  pivot_longer(3:6, values_to = "pearson", names_to = "method") %>%
  mutate(increase = pearson) %>%
  mutate(increase = ifelse(pearson >= 1, TRUE,
                           ifelse(pearson < 1, FALSE, NA))) %>%
  group_by(method, increase) %>%
  count(increase) %>%
  filter(increase == TRUE)

# out degree centrality
# centrality(data_to_adjacency(day_1[1])[[1]] %>% as.matrix, measure = "degree")

out_degree_intersect <- function(adj_1, adj_2) {
  #' Takes two adjacency matrices "adj_1" and "adj_2" and computes their out-degree
  #' centrality for the intersecting symptoms in both adjacency matrices. The function
  #' returns a logical value of whether the symptoms with the highest out-degree centrality
  #' is the same for both adjacency matrices.

  intersection <- intersect(rownames(adj_1), rownames(adj_2))

  out_deg_1 <- adj_1 %>%
    as.matrix() %>%
    graph_from_adjacency_matrix(mode = "directed", weighted = TRUE) %>%
    strength(mode = "out") %>%
    data.frame(symptoms = rownames(adj_1)) %>%
    rename("out_degree" = ".") %>%
    arrange(out_degree) %>%
    filter(symptoms %in% intersection) %>%
    tail(1)

  out_deg_2 <- adj_2 %>%
    as.matrix() %>%
    graph_from_adjacency_matrix(mode = "directed", weighted = TRUE) %>%
    strength(mode = "out") %>%
    data.frame(symptoms = rownames(adj_2)) %>%
    rename("out_degree" = ".") %>%
    arrange(out_degree) %>%
    filter(symptoms %in% intersection) %>%
    tail(1)

  return(out_deg_1$symptoms == out_deg_2$symptoms)
}

compare_out_degree <- function(f,...) {
  # computes the ratio of number of patients having the same symptom with the
  # maximum out-degree centrality, when applying the function f. Could be either
  # data_to_adjacency or one of the pruning methods.
  n_equal <- 0
  for (i in 1:length(day_1)) {
    highest_equal <- out_degree_intersect(f(data_to_adjacency(day_1[i]))[[1]],
                                          f(data_to_adjacency(day_2[i]))[[1]]
    )
    if (highest_equal == TRUE) {
      n_equal <- n_equal + 1
    }

  }

  return(n_equal / length(day_1))
}

compare_out_degree(function(x){x = x})

compare_out_degree <- function(f,...) {
  # computes the ratio of number of patients having the same symptom with the
  # maximum out-degree centrality, when applying the function f. Could be either
  # data_to_adjacency or one of the pruning methods.
  n_equal <- 0
  for (i in 1:length(day_1)) {
    highest_equal <- out_degree_intersect(f(data_to_adjacency(day_1[i]))[[1]] %>%
                                            as_adjacency_matrix(attr = "weight", sparse = FALSE),
                                          f(data_to_adjacency(day_2[i]))[[1]] %>%
                                            as_adjacency_matrix(attr = "weight", sparse = FALSE)
    )
    if (highest_equal == TRUE) {
      n_equal <- n_equal + 1
    }

  }

  return(n_equal / length(day_1))
}

# compare_out_degree(function(x){x = x})
compare_out_degree(pruning_edge_betweenness)
compare_out_degree(pruning_page_rank)
compare_out_degree(pruning_updating_page_rank)
compare_out_degree(pruning_connectivity_kept)


# ------------------------------------------------------------------------------
# Adding noise
# ------------------------------------------------------------------------------

edge_weights_small <- rnorm(5, mean = 70, sd = 10) %>% round(0)
small_adj <- rbind(c(0,edge_weights_small[1],0,0,0),
                   c(edge_weights_small[2],0,edge_weights_small[3],edge_weights_small[4],edge_weights_small[5]),
                   c(0,0,0,0,0),
                   c(0,0,0,0,0),
                   c(0,0,0,0,0)) %>%
  as.matrix()
mean(small_adj[small_adj > 0]) # mean edge weight
small_graph <- graph_from_adjacency_matrix(small_adj, mode = "directed", weighted = TRUE)
V(small_graph)$name <- V(small_graph) %>% as_ids()
small_node_weights <- data_frame(node = V(small_graph)$name %>% as.character(), "node_weight" = rnorm(5, mean = 70, sd = 10) %>% round(0))
mean(small_node_weights$node_weight) # mean node weight

### large graph ###
edge_weights_large <- rnorm(12, mean = 70, sd = 10) %>% round(0)
large_adj <- rbind(c(0,edge_weights_large[1],edge_weights_large[2],edge_weights_large[3],0,
                     edge_weights_large[5],edge_weights_large[6],0,edge_weights_large[8]),
                   c(edge_weights_large[9],0,0,0,0,0.0,0,0,0),
                   c(0,0,0,0,0,0.0,0,0,0),
                   c(0,0,0,0,0,0.0,0,0,0),
                   c(edge_weights_large[4],0,0,0,0,0.0,0,0,0),
                   c(edge_weights_large[12],0,0,0,0,0.0,0,0,0),
                   c(0,0,0,0,0,0.0,0,0,0),
                   c(edge_weights_large[7],0,0,0,0,0.0,0,0,0),
                   c(0,0,0,0,0,0.0,0,edge_weights_large[10],0))

mean(large_adj[large_adj > 0]) # mean edge weight
large_graph <- graph_from_adjacency_matrix(large_adj, mode = "directed", weighted = TRUE)
V(large_graph)$name <- V(large_graph) %>% as_ids()
large_node_weights <- data_frame(node = V(large_graph)$name %>% as.character(), "node_weight" = rnorm(9, mean = 70, sd = 10) %>% round(0))
mean(large_node_weights$node_weight) # mean node weight

amount_of_noise <- function(g, frac_nodes, frac_existing_edges, frac_new_edges){
  #' samples the noise so we can add and remove it and see effects without
  #' resampling new noise for the four noise types. Also samples matrix and
  #' vector indices of where to add the noise to only do this once per repeated
  #' trial. The graph is given by "g", the fractions of nodes, existing edges, and
  #' new edges to add noise to are given by "frac_nodes", "frac_existing_edges",
  #' and "frac_new_edges", respectively.

  adj <- g %>%
    as_adjacency_matrix(attr = "weight", sparse = FALSE)
  diag(adj) <- NA # to make sure we do not add loops

  # finding number of total possible nodes, existing edges, and new edges to add noise to
  n_existing_edges <- length(E(g))
  n_new_edges <- nrow(adj)^2 - length(E(g)) - length(diag(adj)) # no. possible edges to add
  n_nodes <- length(V(g))

  #### adding noise to networks ###
  # finding number of nodes, existing edges, and new edges to add noise to using the given fraction
  n_noisy_nodes <- ceiling(n_nodes * frac_nodes)
  n_noisy_edges_existing <- ceiling(n_existing_edges * frac_existing_edges)
  n_noisy_edges_new <- ceiling(n_new_edges * frac_new_edges)

  # sample of indices to add noise to
  noise_inds_nodes <- sample(1:n_nodes, n_noisy_nodes) # indices of the nodes to add noise to
  noise_inds_exist <- sample(which(t(adj) > 0), n_noisy_edges_existing) # indices of the existing edges to add noise to
  noise_inds_new <- sample(which(t(adj) == 0), n_noisy_edges_new) # indices of the potential edges to add as noise
  diag(adj) <- 0 # add back 0 on diagonal when we have chosen where to add new edges

  epsilon_V <- rnorm(n_noisy_nodes)
  epsilon_E <- rnorm(n_noisy_edges_existing)
  epsilon_E_hat <- rnorm(n_noisy_edges_new)

  return(list(epsilon_V, epsilon_E, epsilon_E_hat,
              noise_inds_nodes, noise_inds_exist, noise_inds_new))
}

adding_noise <- function(g, node_weights, file_path = "",
                         frac_nodes, frac_existing_edges, frac_new_edges,
                         lambda_V, lambda_E, lambda_E_hat,
                         epsilon_V, epsilon_E, epsilon_E_hat,
                         noise_inds_nodes, noise_inds_exist, noise_inds_new){
  #' adds noise to nodes and existing edges, and adds new edges to the graph
  #' "g" given the original "node_weights" and the "file_path" of the
  #' patient (used for renaming a data frame). The fractions of nodes, existing edges, and
  #' new edges to add noise to are given by "frac_nodes", "frac_existing_edges",
  #' and "frac_new_edges", respectively. The amounts of noise for these are given
  #' by "noise_inds_nodes", "noise_inds_exist", and "noise_inds_new", respectively, and
  #' the remaining arguments are according to the notation in the thesis.

  # setting up adjacency matrix and node weights data frame
  adj <- g %>%
    as_adjacency_matrix(attr = "weight", sparse = FALSE)
  no_change_adj <- adj

  df_node_w <- node_weights[[2]] %>%
    data_frame() %>%
    mutate(node = V(g)$name %>% as.character()) %>%
    rename("node_weights" = ".") %>%
    relocate(node)

  # finding number of total possible nodes, existing edges, and new edges to add noise to
  n_existing_edges <- length(E(g))
  n_new_edges <- nrow(adj)^2 - length(E(g)) - length(diag(adj)) # no. possible edges to add
  n_nodes <- length(V(g))

  #### adding noise to networks ###

  # column id is the remainder when dividing with # nodes (we only find the indices
  ## in the list and not the actual matrix from the above)
  col_id_exist <- noise_inds_exist %% ncol(adj)
  col_id_exist[col_id_exist == 0] <- ncol(adj) # accounts for remainder 0 if in last col

  col_id_new <- noise_inds_new %% ncol(adj)
  col_id_new[col_id_new == 0] <- ncol(adj) # accounts for remainder 0 if in last col

  ## row id is 1 + the quotient when subtracting the remainder and dividing with # nodes
  row_id_exist <- (noise_inds_exist - col_id_exist) / ncol(adj) + 1
  row_id_new <- (noise_inds_new - col_id_new) / ncol(adj) + 1

  ## final ids
  ids_exist <- cbind(row_id_exist, col_id_exist)
  ids_new <- cbind(row_id_new, col_id_new)

  ## adding noise to node weights
  noisy_node_weights <- df_node_w # before we add noise
  v_noisy <- df_node_w$node_weights[noise_inds_nodes] * (1 + lambda_V * epsilon_V)
  index_neg <- noise_inds_nodes[which(v_noisy <= 0)]
  # if new noisy weights are negative or zero, don't add noise
  if (length(index_neg) != 0) {

    v_noisy[v_noisy <= 0] <- df_node_w$node_weights[index_neg]
  }

  # if new noisy weights are greater than 100, cutoff at 100
  v_noisy[v_noisy > 100] <- 100
  noisy_node_weights$node_weights[noise_inds_nodes] <- v_noisy # replacing old weights with new noisy weights
  absolute_node_noise <- abs(df_node_w$node_weights - noisy_node_weights$node_weights)

  ## adding noise to existing edge weights
  noisy_adj <- adj # before we add noise
  w_noisy <- adj[ids_exist] * (1 + lambda_E * epsilon_E)
  row_neg <- ids_exist[which(w_noisy <= 0), 1]
  col_neg <- ids_exist[which(w_noisy <= 0), 2]
  indices_neg <- data_frame(row_neg, col_neg)
  # if new noisy weights are negative or zero, don't add noise
  w_noisy[w_noisy <= 0] <- adj[cbind(indices_neg$row_neg, indices_neg$col_neg)]
  # if new noisy weights are greater than 100, cutoff at 100
  w_noisy[w_noisy > 100] <- 100
  noisy_adj[ids_exist] <- w_noisy # replacing old weights with new noisy weights
  absolute_edge_noise <- abs(adj - noisy_adj)

  ## adding noise in the form of new edge weights
  noisy_adj_before_new <- noisy_adj # before we add new noisy edges to the noisy adj matrix
  w_noisy_new <- lambda_E_hat * mean(no_change_adj[no_change_adj > 0]) * (1 + lambda_E * epsilon_E_hat)
  row_neg <- ids_new[which(w_noisy_new <= 0), 1]
  col_neg <- ids_new[which(w_noisy_new <= 0), 2]
  indices_neg <- data_frame(row_neg, col_neg)
  # if new noisy weights are negative or zero, don't add noise
  w_noisy_new[w_noisy <= 0] <- adj[cbind(indices_neg$row_neg, indices_neg$col_neg)]
  # if new noisy weights are greater than 100, cutoff at 100
  w_noisy_new[w_noisy > 100] <- 100
  noisy_adj[ids_new] <- w_noisy_new # replacing old edge weights w new edge weights for the selected new edges
  absolute_new_edge_noise <- abs(w_noisy_new)

  return(list(noisy_adj, noisy_node_weights, no_change_adj, df_node_w,
              absolute_node_noise, absolute_edge_noise, absolute_new_edge_noise))
}

# Plotting original small and large synthetic networks

plot_weighted_network(small_graph, "",
                      small_node_weights$node_weight * 30 / 100,
                      E(small_graph)$weight * 30 / 100,
                      text_size = 4)

plot_weighted_network(large_graph, "",
                      large_node_weights$node_weight * 30 / 100,
                      E(large_graph)$weight * 30 / 100,
                      text_size = 4)
# all noise
noisy_small_all <- adding_noise(small_graph, small_node_weights, "", 5/5,5/5,5/15,0.2,0.2,0.2, ,
                                epsilon_V, epsilon_E, epsilon_E_hat)
plot_weighted_network(noisy_small_all[[1]] %>% graph_from_adjacency_matrix(mode = "directed", weighted = TRUE),
                      "",
                      noisy_small_all[[2]]$node_weights / 100 * 30,
                      E(noisy_small_all[[1]] %>% graph_from_adjacency_matrix(mode = "directed",
                                                                             weighted = TRUE))$weight,
                      text_size = 4)

## ------------------------------ ##
## Jaccard coefficients for noise ##
## ------------------------------ ##

jaccard_tables_noise <- function(adj_node_w_noisy, adj_node_w_og) {
  #' Computes a table with the Jaccard index computed for the graphs represented by
  #' "adj_node_w_noisy" and "adj_node_w_og" for all four methods and the original
  #' synthetic network.

  #### prune noisy graph ####
  # naive edge betweenness
  pruning_edge_betweenness_old <- pruning_edge_betweenness(adj_node_w_noisy)[[1]] %>% E() %>% as_ids()

  # edge betweenness with PageRank
  pruning_page_rank_old <- pruning_page_rank(adj_node_w_noisy)[[1]] %>% E() %>% as_ids()

  # Updated PageRank
  pruning_updating_page_rank_old <- pruning_updating_page_rank(adj_node_w_noisy)[[1]] %>% E() %>% as_ids()

  # Brute force
  pruning_connectivity_kept_old <- pruning_connectivity_kept(adj_node_w_noisy)[[1]] %>% E() %>% as_ids()

  #### Jaccard ####
  jaccard <- list("original" = jaccard(graph_from_adjacency_matrix(adj_node_w_noisy[[1]] %>% as.matrix(),
                                                                   mode = "directed",
                                                                   weighted = TRUE) %>% E() %>% as_ids(),
                                       graph_from_adjacency_matrix(adj_node_w_og[[1]] %>% as.matrix(),
                                                                   mode = "directed",
                                                                   weighted = TRUE) %>% E() %>% as_ids()),
                  "edge_b" = jaccard(pruning_edge_betweenness_old,
                                     graph_from_adjacency_matrix(adj_node_w_og[[1]] %>%
                                                                   as.matrix(),
                                                                 mode = "directed",
                                                                 weighted = TRUE) %>% E() %>% as_ids()),
                  "page_rank" = jaccard(pruning_page_rank_old,
                                        graph_from_adjacency_matrix(adj_node_w_og[[1]] %>% as.matrix(),
                                                                    mode = "directed",
                                                                    weighted = TRUE) %>% E() %>% as_ids()),
                  "updated_pr" = jaccard(pruning_updating_page_rank_old,
                                         graph_from_adjacency_matrix(adj_node_w_og[[1]] %>% as.matrix(),
                                                                     mode = "directed",
                                                                     weighted = TRUE) %>% E() %>% as_ids()),
                  "connectivity_kept" = jaccard(pruning_connectivity_kept_old,
                                          graph_from_adjacency_matrix(adj_node_w_og[[1]] %>% as.matrix(),
                                                                      mode = "directed",
                                                                      weighted = TRUE) %>% E() %>% as_ids())
  )

  #### creating a table ####
  library(data.table)
  jaccard_to_table <- map(jaccard, as.data.table)
  jaccard_table <- rbindlist(jaccard_to_table, fill = TRUE, idcol = TRUE) %>%
    rename("method" = ".id", "jaccard_g" = "V1") %>%
    mutate(jaccard_g = jaccard_g %>% signif(3))

  return(jaccard_table %>% data_frame())
}

## ---------------------------------- ##
## Correlation coefficients for noise ##
## ---------------------------------- ##

joining_tables <- function(graph, node_weights, network = "",
                           frac_nodes, frac_existing_edges, frac_new_edges,
                           lambda_v, lambda_E, lambda_E_hat,
                           epsilon_V, epsilon_E, epsilon_E_hat,
                           noise_inds_nodes, noise_inds_exist, noise_inds_new) {
  #' creates a table of the Pearson, Spearman, and Jaccard similarity coefficients given the
  #' "graph", "node_weights", file path "network", "frac_nodes", "frac_existing_edges",
  #' and "frac_new_edges" to add noise to, and the level of noise "lambda_v", "lambda_E",
  #' "lambda_E_hat" for the nodes, existing edges, and new edges, respectively.

  tables_list <- list()
  for (i in 1:5) {
    # subset noise and indices to choose same random noise for each of the four combinations to easier
    # see the differences/effects of adding noise like that
    noise_res <- adding_noise(graph, node_weights, network,
                              frac_nodes, frac_existing_edges, frac_new_edges,
                              lambda_v, lambda_E, lambda_E_hat,
                              epsilon_V[[i]], epsilon_E[[i]], epsilon_E_hat[[i]],
                              noise_inds_nodes[[i]], noise_inds_exist[[i]], noise_inds_new[[i]])

    assign(paste0("noise_adj", i), noise_res[[1]])
    assign(paste0("noise_node", i), noise_res[[2]])
    assign(paste0("abs_node", i), noise_res[[5]])
    assign(paste0("abs_edge", i), noise_res[[6]])
    assign(paste0("abs_new_edge", i), noise_res[[7]])

    assign(paste0("df_cor", i), correlation_tables_intersect(list(noise_res[[1]], noise_res[[2]]),
                                                             list(noise_res[[3]], noise_res[[4]]),
                                                             method = "pearson",
                                                             correlation))
    assign(paste0("df_cor_ranked", i), correlation_tables_intersect(list(noise_res[[1]], noise_res[[2]]),
                                                                    list(noise_res[[3]], noise_res[[4]]),
                                                                    method = "pearson",
                                                                    correlation_ranked))

    assign(paste0("df_j", i), jaccard_tables_noise(list(noise_res[[1]], noise_res[[2]]),
                                                   list(noise_res[[3]], noise_res[[4]])))
  }

  tables_list_noise_adj <- list(noise_adj1, noise_adj2, noise_adj3, noise_adj4, noise_adj5)
  tables_list_noise_node <- list(noise_node1, noise_node2, noise_node3, noise_node4, noise_node5)

  tables_list_abs_node <- list(abs_node1, abs_node2, abs_node3, abs_node4, abs_node5)
  tables_list_abs_edge <- list(abs_edge1, abs_edge2, abs_edge3, abs_edge4, abs_edge5)
  tables_list_abs_new_edge <- list(abs_new_edge1, abs_new_edge2, abs_new_edge3, abs_new_edge4, abs_new_edge5)

  tables_list_cor <- list(df_cor1,df_cor2,df_cor3,df_cor4,df_cor5)
  tables_list_cor_ranked <- list(df_cor_ranked1,df_cor_ranked2,df_cor_ranked3,df_cor_ranked4,df_cor_ranked5)
  tables_list_j <- list(df_j1,df_j2,df_j3,df_j4,df_j5)

  # Creating data frame of all similarity indices
  table_cor <- plyr::join_all(tables_list_cor, by = "method")
  table_cor_ranked <- plyr::join_all(tables_list_cor_ranked, by = "method")
  table_j <- plyr::join_all(tables_list_j, by = "method")

  # Changing column names to be able to find row means
  colnames <- names(table_cor)
  for (i in 2:ncol(table_cor)) {
    colnames[i] <- paste0(colnames[i], "_", i - 1)
  }
  colnames(table_cor) <- colnames

  # Changing column names to be able to find row means
  colnames <- names(table_cor_ranked)
  for (i in 2:ncol(table_cor_ranked)) {
    colnames[i] <- paste0(colnames[i], "_", i - 1)
  }
  colnames(table_cor_ranked) <- colnames

  # Changing column names to be able to find row means
  colnames <- names(table_j)
  for (i in 2:ncol(table_j)) {
    colnames[i] <- paste0(colnames[i], "_", i - 1)
  }
  colnames(table_j) <- colnames

  return(list(table_cor, table_cor_ranked, table_j,
              tables_list_abs_node, tables_list_abs_edge, tables_list_abs_new_edge,
              tables_list_noise_adj, tables_list_noise_node))
}

noise_results_table <- function(results) {
  #' joins the Pearson, Spearman, and Jaccard results in a single table, given the noise
  #' results in the argument "results".

  pearson <- results[[1]] %>%
    select(method, starts_with("cor_g")) %>%
    pivot_longer(2:6,values_to = "cor", names_to = "person") %>%
    pivot_wider(names_from = "method", values_from = "cor") %>%
    pivot_longer(2:6, values_to = "pearson", names_to = "method") %>%
    group_by(method) %>%
    summarize(mean(pearson)) %>%
    slice(match(c("original", "edge_b", "page_rank", "updated_pr", "connectivity_kept"), method))

  spearman <- results[[2]] %>%
    select(method, starts_with("cor_g")) %>%
    pivot_longer(2:6,values_to = "cor", names_to = "person") %>%
    pivot_wider(names_from = "method", values_from = "cor") %>%
    pivot_longer(2:6, values_to = "pearson", names_to = "method") %>%
    group_by(method) %>%
    summarize(mean(pearson)) %>%
    slice(match(c("original", "edge_b", "page_rank", "updated_pr", "connectivity_kept"), method))

  jaccard <- results[[3]] %>%
    select(method, starts_with("jaccard_g")) %>%
    pivot_longer(2:6,values_to = "jaccard", names_to = "person") %>%
    pivot_wider(names_from = "method", values_from = "jaccard") %>%
    pivot_longer(2:6, values_to = "jaccard", names_to = "method") %>%
    group_by(method) %>%
    summarize(mean(jaccard)) %>%
    slice(match(c("original", "edge_b", "page_rank", "updated_pr", "connectivity_kept"), method))

  pearson %>%
    inner_join(spearman, by = "method") %>%
    inner_join(jaccard, by = "method") %>%
    return()
}

abs_mean_noise <- function(results){
  #' computes the absolute mean noise added to nodes, existing edges, and new
  #' edges, respectively, and returns the results in a list in that order.

  # abs node noise
  nodes <- map(results[[4]], mean) %>% map(digits = 2, round)
  # abs existing edge noise, need to compute mean ourselves bc 0s
  existing_edges <- map(results[[5]], sum) %>% map(~. / 5) %>% map(digits = 2, round)
  # abs weight for new edges
  new_edges <- map(results[[6]], mean) %>% map(digits = 2, round)

  list(nodes, existing_edges, new_edges) %>%
    return()
}

get_abs_means_5 <- function(results) {
  #' computes the means of the mean absolute noise values from the 5 repeated noise trials

  nodes <- abs_mean_noise(results)[[1]] %>% as.numeric() %>% mean() %>% round(2)
  existing_edges <- abs_mean_noise(results)[[2]] %>% as.numeric() %>% mean() %>% round(2)
  new_edge_weights <- abs_mean_noise(results)[[3]] %>% as.numeric() %>% mean() %>% round(2)

  return(c(nodes, existing_edges, new_edge_weights))
}

plot_noisy_network <- function(results) {
  #' plots the 5 noisy networks given the results
  for(i in 1:5){
    g <- results[[7]][[i]] %>%
      graph_from_adjacency_matrix(mode = "directed", weighted = TRUE)
    g %>%
      plot_weighted_network(title = "", results[[8]][[i]]$node_weights * 0.3, E(g)$weight * 0.3, text_size = 4) %>%
      print()
  }
}
plot_noisy_network_unweighted <- function(results) {
  #' plots the 5 noisy networks given the results
  for(i in 1:5){
    g <- results[[7]][[i]] %>%
      graph_from_adjacency_matrix(mode = "directed", weighted = TRUE)
    g %>%
      plot_weighted_network(title = "", 30, 30, text_size = 4) %>%
      print()
  }
}

#' computing the Pearson, Spearman, and Jaccard coefficients for 5 noisy networks of each of the
#' synthetic networks below, with varying noise
results_small <- joining_tables(small_graph, small_node_weights, "",
                                frac_nodes = 5/5, frac_existing_edges = 5/5, frac_new_edges = 7/15,
                                lambda_v = 0.1, lambda_E = 0.1, lambda_E_hat = 0.1)
# results_small <- joining_tables(small_graph, small_node_weights, "",
#                                 frac_nodes = 0/5, frac_existing_edges = 0/5, frac_new_edges = 10/15,
#                                 lambda_v = 0.1, lambda_E = 0.1, lambda_E_hat = 0.1)

sampling_noise_indices_5 <- function(g, frac_nodes, frac_existing_edges, frac_new_edges) {
  #' samples the noise and indices to add noise to given a graph "g" and the fractions of
  #' nodes, existing edges, and new edges to add the noise to according to the arguments
  #' "frac_nodes", "frac_existing_edges", and "frac_new_edges".

  # sample 5 times
  noise_amount_large_1 <- amount_of_noise(g, frac_nodes, frac_existing_edges, frac_new_edges)
  noise_amount_large_2 <- amount_of_noise(g, frac_nodes, frac_existing_edges, frac_new_edges)
  noise_amount_large_3 <- amount_of_noise(g, frac_nodes, frac_existing_edges, frac_new_edges)
  noise_amount_large_4 <- amount_of_noise(g, frac_nodes, frac_existing_edges, frac_new_edges)
  noise_amount_large_5 <- amount_of_noise(g, frac_nodes, frac_existing_edges, frac_new_edges)

  # collect in list to return it
  noise_amount_large <- list(list(noise_amount_large_1[[1]], noise_amount_large_2[[1]], noise_amount_large_3[[1]],
                                  noise_amount_large_4[[1]], noise_amount_large_5[[1]]),
                             list(noise_amount_large_1[[2]], noise_amount_large_2[[2]], noise_amount_large_3[[2]],
                                  noise_amount_large_4[[2]], noise_amount_large_5[[2]]),
                             list(noise_amount_large_1[[3]], noise_amount_large_2[[3]], noise_amount_large_3[[3]],
                                  noise_amount_large_4[[3]], noise_amount_large_5[[3]]),
                             list(noise_amount_large_1[[4]], noise_amount_large_2[[4]], noise_amount_large_3[[4]],
                                  noise_amount_large_4[[4]], noise_amount_large_5[[4]]),
                             list(noise_amount_large_1[[5]], noise_amount_large_2[[5]], noise_amount_large_3[[5]],
                                  noise_amount_large_4[[5]], noise_amount_large_5[[5]]),
                             list(noise_amount_large_1[[6]], noise_amount_large_2[[6]], noise_amount_large_3[[6]],
                                  noise_amount_large_4[[6]], noise_amount_large_5[[6]])
  )
  return(noise_amount_large)
}

### Examples of how noise was added + plots ###
noise_amount_large <- sampling_noise_indices_5(large_graph,
                                               frac_nodes = 9/9, frac_existing_edges = 11/11, frac_new_edges = 39/61)

results_large <- joining_tables(large_graph, large_node_weights, "",
                                frac_nodes = 9/9, frac_existing_edges = 11/11, frac_new_edges = 39/61,
                                lambda_v = 0.2, lambda_E = 0.2, lambda_E_hat = 0.9,
                                epsilon_V = noise_amount_large[[1]],
                                epsilon_E = noise_amount_large[[2]],
                                epsilon_E_hat = noise_amount_large[[3]],
                                noise_inds_nodes = noise_amount_large[[4]],
                                # noise_inds_nodes = rep(NULL,5),
                                noise_inds_exist = noise_amount_large[[5]],
                                # noise_inds_exist = rep(NULL,5),
                                noise_inds_new = noise_amount_large[[6]]
)

noise_amount_small <- sampling_noise_indices_5(small_graph,
                                               frac_nodes = 5/5, frac_existing_edges = 5/5, frac_new_edges = 10/15)

results_small <- joining_tables(small_graph, small_node_weights, "",
                                frac_nodes = 5/5, frac_existing_edges = 5/5, frac_new_edges = 10/15,
                                lambda_v = 0.1, lambda_E = 0.1, lambda_E_hat = 0.1,
                                epsilon_V = noise_amount_small[[1]],
                                epsilon_E = noise_amount_small[[2]],
                                epsilon_E_hat = noise_amount_small[[3]],
                                noise_inds_nodes = noise_amount_small[[4]],
                                # noise_inds_nodes = rep(NULL,5),
                                noise_inds_exist = noise_amount_small[[5]],
                                # noise_inds_exist = rep(NULL,5),
                                noise_inds_new = noise_amount_small[[6]]
)

# Pearson, Spearman, Jaccard for original and pruned networks
noise_results_table(results_large)
noise_results_table(results_small)
noise_results_table(results_large_1)

# Absolute mean node noise, existing edge noise, and weights of new edges added as noise
get_abs_means_5(results_small)
get_abs_means_5(results_large)

# Check out-degree
map(results_large[[7]], adj_2 = large_adj,out_degree_noise)

# Plots of 1 of the 5 noisy networks (examples)
plot_noisy_network(results_large)
plot_noisy_network_unweighted(results_small)

plot_weighted_network(large_graph, "",
                      large_node_weights$node_weight * 30 / 100,
                      E(large_graph)$weight * 30 / 100,
                      text_size = 4)
