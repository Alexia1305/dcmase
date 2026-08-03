library(Matrix)
library(aricode)
library(reticulate)
source("R/comdetmethods.R")
library(igraph)
library(mclust)

####################### Load DATA ##########################################################

build_AUCS <- function(edge_file, node_file){

  ############################
  # Lecture des données
  ############################
  
  edges <- read.csv(
    edge_file,
    header = TRUE,
    stringsAsFactors = FALSE
  )
  
  nodes <- read.csv(
    node_file,
    header = TRUE,
    stringsAsFactors = FALSE
  )


  ############################
  # Liste des noeuds
  ############################
  
  node_names <- nodes$node
  N <- length(node_names)
  
  node_index <- setNames(1:N, node_names)


  ############################
  # Matrices d'adjacence
  ############################
  
  layers <- unique(edges$layer)
  
  adjacency_list <- list()


  for(layer_name in layers){
    
    A <- Matrix(
      0,
      nrow = N,
      ncol = N,
      sparse = TRUE
    )
    
    rownames(A) <- node_names
    colnames(A) <- node_names
    
    
    e <- edges[edges$layer == layer_name, ]
    
    
    for(i in 1:nrow(e)){
      
      s <- e$source[i]
      t <- e$target[i]
      
      A[node_index[s], node_index[t]] <- 1
      A[node_index[t], node_index[s]] <- 1
    }
    
    
    adjacency_list[[layer_name]] <- A
  }
  names(adjacency_list) <- NULL

  ############################
  # Groupes des noeuds
  ############################
  
  groups <- as.numeric(sub("G", "", nodes$group))
 
  ############################
  # Retour
  ############################
  
  return(
    list(
      A = adjacency_list,
      labels = groups
    )
  )
}

build_cora_multilayer <- function(content_path, cites_path, k = 20) {
  
  # -------------------------
  # 1. Load content file
  # -------------------------
  content <- read.table(content_path,
                        header = FALSE,
                        stringsAsFactors = FALSE)
  
  paper_id <- content[, 1]
  labels_raw <- content[, ncol(content)]
  
  X <- as.matrix(content[, 2:(ncol(content) - 1)])
  rownames(X) <- paper_id
  
  # -------------------------
  # KEEP ONLY 3 CLASSES
  # -------------------------
  keep_classes <- c(
    "Genetic_Algorithms",
    "Neural_Networks",
    "Probabilistic_Methods"
  )
  
  keep_idx <- which(labels_raw %in% keep_classes)
  
  paper_id <- paper_id[keep_idx]
  labels_raw <- labels_raw[keep_idx]
  X <- X[keep_idx, ]
  
  # convert labels to 1..3
  labels <- match(labels_raw, keep_classes)
  
  n <- length(paper_id)
  
  # -------------------------
  # 2. CITATION LAYER
  # -------------------------
  cites <- read.table(cites_path,
                      header = FALSE,
                      stringsAsFactors = FALSE)
  
  colnames(cites) <- c("cited", "citing")
  
  # filter edges to kept nodes
paper_id <- as.character(paper_id)
cites$citing <- as.character(cites$citing)
cites$cited <- as.character(cites$cited)

adj_citation <- matrix(0, n, n,
                       dimnames = list(paper_id, paper_id))

valid <- cites$citing %in% paper_id &
         cites$cited %in% paper_id

cites <- cites[valid, ]

adj_citation[cbind(cites$citing, cites$cited)] <- 1
adj_citation[cbind(cites$cited, cites$citing)] <- 1
  # -------------------------
  # 3. SIMILARITY LAYER
  # -------------------------
  norm_X <- sqrt(rowSums(X^2))
  norm_X[norm_X == 0] <- 1  # avoid division by zero
  
  X_norm <- X / norm_X
  
  sim <- X_norm %*% t(X_norm)
  diag(sim) <- 0
  
  # -------------------------
  # 4. kNN GRAPH
  # -------------------------
  adj_similarity <- matrix(0, n, n)
  rownames(adj_similarity) <- paper_id
  colnames(adj_similarity) <- paper_id
  
  for (i in 1:n) {
    top_k <- order(sim[i, ], decreasing = TRUE)[1:min(k, n)]
    adj_similarity[i, top_k] <- 1
  }
  
  # symmetrize
 adj_similarity <- pmax( adj_similarity, t( adj_similarity))
  
  # -------------------------
  # 5. RETURN
  # -------------------------
    A <- list(
  adj_citation,
  adj_similarity)
  return(list(
    A = A ,
    labels           = labels
  ))
}

build_citeseer_multilayer <- function(content_path, cites_path, k = 20) {
  
  # -------------------------
  # 1. Load content file
  # -------------------------
  content <- read.table(content_path,
                        header = FALSE,
                        stringsAsFactors = FALSE)
  
  paper_id <- content[, 1]
  labels_raw <- content[, ncol(content)]
  
  X <- as.matrix(content[, 2:(ncol(content) - 1)])
  rownames(X) <- paper_id
  
  classes <- c(
    "Agents",
			"AI",
			"DB",
			"IR",
			"ML",
			"HCI"
  )
  # convert labels to 1..3
  labels <- match(labels_raw, classes)
  
  n <- length(paper_id)
  
  # -------------------------
  # 2. CITATION LAYER
  # -------------------------
  cites <- read.table(cites_path,
                      header = FALSE,
                      stringsAsFactors = FALSE)
  
  colnames(cites) <- c("cited", "citing")
  
  # filter edges to kept nodes
paper_id <- as.character(paper_id)
cites$citing <- as.character(cites$citing)
cites$cited <- as.character(cites$cited)

adj_citation <- matrix(0, n, n,
                       dimnames = list(paper_id, paper_id))

valid <- cites$citing %in% paper_id &
         cites$cited %in% paper_id

cites <- cites[valid, ]

adj_citation[cbind(cites$citing, cites$cited)] <- 1
adj_citation[cbind(cites$cited, cites$citing)] <- 1
  # -------------------------
  # 3. SIMILARITY LAYER
  # -------------------------
  norm_X <- sqrt(rowSums(X^2))
  norm_X[norm_X == 0] <- 1  # avoid division by zero
  
  X_norm <- X / norm_X
  
  sim <- X_norm %*% t(X_norm)
  diag(sim) <- 0
  
  # -------------------------
  # 4. kNN GRAPH
  # -------------------------
  adj_similarity <- matrix(0, n, n)
  rownames(adj_similarity) <- paper_id
  colnames(adj_similarity) <- paper_id
  
  for (i in 1:n) {
    top_k <- order(sim[i, ], decreasing = TRUE)[1:min(k, n)]
    adj_similarity[i, top_k] <- 1
  }
  
  # symmetrize
   adj_similarity <- pmax( adj_similarity, t( adj_similarity))
  
  # -------------------------
  # 5. RETURN
  # -------------------------

  A <- list(
  adj_citation,
  adj_similarity)
  return(list(
    A = A ,
    labels           = labels
  ))
}

build_UCI<- function(path, k = 20) {
  
  files <- c("mfeat-fou", "mfeat-fac", "mfeat-kar",
             "mfeat-pix", "mfeat-zer", "mfeat-mor")
  
  # labels: 200 instances per class (0–9)
  true_labels <- rep(0:9, each = 200)
  
  adjacency_list <- lapply(files, function(f) {
    
    file_path <- file.path(path, f)
    
    # load data
    X <- as.matrix(read.table(file_path))
    
    n <- nrow(X)
    
    # Euclidean distances
    dist_mat <- as.matrix(dist(X, method = "euclidean"))
    
    # k-NN adjacency (directed first)
    A <- matrix(0, n, n)
    
    for (i in 1:n) {
      nn <- order(dist_mat[i, ])[2:(k + 1)]
      A[i, nn] <- 1
    }
    
    # make graph non-oriented (symmetrize)
   
    A <- pmax(A, t(A))
    return(A)
  })
  
  
  return(list(
    A = adjacency_list,
    labels = true_labels
  ))
}

build_CBCL<- function(edges_path, labels_path) {
  
  # --- Lecture des labels ---
  labels_df <- read.table(labels_path, header = FALSE, 
                           col.names = c("nodeID", "label"))
  labels_df <- labels_df[order(labels_df$nodeID), ]  # s'assurer de l'ordre
  labels <- labels_df$label
  
  n_nodes <- max(labels_df$nodeID)
  
  
  
  edges_df <- read.table(edges_path, header = FALSE,
                          col.names = c("layer", "i", "j"))
  edges_df$weight <- 1
  
  
  layer_ids <- sort(unique(edges_df$layer))
  
  # --- Construction d'une matrice d'adjacence par couche ---
  A_list <- lapply(layer_ids, function(l) {
    A <- matrix(0, nrow = n_nodes, ncol = n_nodes)
    edges_l <- edges_df[edges_df$layer == l, ]
    
    for (k in seq_len(nrow(edges_l))) {
      i <- edges_l$i[k]
      j <- edges_l$j[k]
      A[i, j] <- 1
      A[j, i] <- 1  
    }
    
    A
  })
  
 
  
   
  return(list(
    A = A_list,
    labels = labels
  ))
}

##################### Methods #############################################################
run_all_methods <- function(Adj_list, truecoms) {
  set.seed(42)
  idx <- !is.na(truecoms)
  K <- length(unique(truecoms[idx]))
  
  # methods_to_run <- c("graph-tool","frost-mf","frost-us","frost-dcmase",
  #                     "us","mf","lmfo","dcmase","ave_spherical",
  #                     "sq-bias-adjusted","mase-spherical")
   methods_to_run <- c("graph-tool","frost-us","dcmase","ave_spherical",
                       "sq-bias-adjusted","mase-spherical","lmfo")
  results <- lapply(methods_to_run, function(method) {
    print(method)
    
    labels <- allmethods(Adj_list, K, method = method)
    
    res<-list(
      NMI = NMI(truecoms[idx], as.vector(labels)[idx]),
      errorRate = classError(labels[idx], truecoms[idx])$errorRate,
      ARI = ARI(truecoms[idx], as.vector(labels)[idx])
    )
    print(res)
  })
  
  # names(results) <- c("graph-tool","FROST_MF","FROST_US","FROST_DCMASE",
  #                     "US","MF","OLMF","DC_MASE","Sum A",
  #                     "S-A^2-Bias-adj","MASE")
  names(results) <- c("graph-tool","frost-us","dcmase","ave_spherical",
                       "sq-bias-adjusted","mase-spherical")
  return(results)
}

library(igraph)

multilayer_properties <- function(adj_list) {
  
  L <- length(adj_list)
  N <- nrow(adj_list[[1]])
  
  cat("Number of nodes:", N, "\n")
  cat("Number of layers:", L, "\n\n")
  
  results <- data.frame(
    layer = 1:L,
    edges = NA,
    density = NA,
    avg_degree = NA,
    degree_sd = NA,
    max_degree = NA,
    isolated_nodes = NA,
    clustering = NA,
    components = NA
  )
  
  degrees <- list()
  
  for (l in 1:L) {
    
    A <- adj_list[[l]]
    
    # Remove diagonal
    diag(A) <- 0
    
    g <- graph_from_adjacency_matrix(
      A,
      mode = "undirected"
    )
    
    deg <- degree(g)
    degrees[[l]] <- deg
    
    results$edges[l] <- ecount(g)
    
    # Density
    results$density[l] <- edge_density(g)
    
    # Average degree
    results$avg_degree[l] <- mean(deg)
    
    # Degree heterogeneity
    results$degree_sd[l] <- sd(deg)
    
    results$max_degree[l] <- max(deg)
    
    # Isolated nodes
    results$isolated_nodes[l] <- sum(deg == 0)/N
    
    # Clustering coefficient
    results$clustering[l] <- transitivity(
      g,
      type = "global"
    )
    
    # Connected components
    results$components[l] <- components(g)$no
  }
  
  
  cat("---- Layer statistics ----\n")
  print(results)
  
  
  cat("\n---- Average properties across layers ----\n")
  
  summary <- data.frame(
    measure = c(
      "Density",
      "Average degree",
      "Degree std",
      "Max degree",
      "Isolated nodes",
      "Clustering",
      "Number of components"
    ),
    mean = c(
      mean(results$density),
      mean(results$avg_degree),
      mean(results$degree_sd),
      mean(results$max_degree),
      mean(results$isolated_nodes),
      mean(results$clustering),
      mean(results$components)
    ),
    sd = c(
      sd(results$density),
      sd(results$avg_degree),
      sd(results$degree_sd),
      sd(results$max_degree),
      sd(results$isolated_nodes),
      sd(results$clustering),
      sd(results$components)
    )
  )
  
  print(summary)
  
  
  # Similarity between layers
  cat("\n---- Layer similarity ----\n")
  
  if(L > 1){
    sim <- matrix(0,L,L)
    
    for(i in 1:L){
      for(j in 1:L){
        
        Ai <- adj_list[[i]]
        Aj <- adj_list[[j]]
        
        sim[i,j] <- sum(Ai == Aj)/(N*N)
      }
    }
    
    print(round(sim,3))
    
    cat("\nAverage layer similarity:",
        mean(sim[upper.tri(sim)]),
        "\n")
  }
  
  
  return(
    list(
      layer_statistics = results,
      summary = summary,
      degrees = degrees
    )
  )
}



######################################## RESULTS ############################

#data <- build_cora_multilayer("Data/cora/cora.content", "Data/cora/cora.cites")
#data <- build_UCI("Data/UCI",k=20)
data <- build_citeseer_multilayer("Data/citeseer/citeseer.content", "Data/citeseer/citeseer.cites")
#data <- build_AUCS("Data/AUCS/aucs_edgelist.txt","Data/AUCS/aucs_nodelist.txt")
#data <- build_CBCL("Data/CBCL/multiplex_edges.txt","Data/CBCL/labels.txt")
properties <- multilayer_properties(data$A)
results<-run_all_methods(data$A[[1]], data$labels)
