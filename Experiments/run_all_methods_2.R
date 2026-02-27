library(Matrix)
library(igraph)
library(mclust)
library(reticulate)
use_python("C:/Users/alexi/anaconda3/python.exe", required = TRUE)
source("R/Codes_Spectral_Matrix_Paul_Chen_AOS_2020.r")
source("R/comdet-dcmase.R")
source("R/comdetmethods.R")
source("R/dcmase.R")
source("R/SpectralMethods.R")
source("R/run_graph_tool.R")
source("R/make_ggplot.R")

sys_path <- "C:/Users/alexi/Documents/MDCBM/Python"
sys <- import("sys")
sys$path <- c(sys$path, sys_path)
frost <- import("frost.frost_multilayer")
np <- import("numpy")


run_simulations <- function(sim_setting, parameters, repetitions = 20) {
  
  results <- lapply(1:repetitions, function(seed) {
    cat("Seed:", seed, "\n")
    generate_data <- sim_setting(parameters, seed)
    run_all_methods(generate_data$Adj_list, generate_data$truecom)
  })
  
  df_res <- data.frame(Reduce(rbind, results))
  rownames(df_res) <- 1:repetitions
  return(df_res)
}

iterate_parameters <- function(sim_setting, parameters_list, param_iter, 
                               repetitions = 20) {
  df_res <- lapply(1:length(parameters_list),  function(i) {
    cat("Running parameter ", param_iter[[i]], "...\n", sep = "")
    sim_res <- run_simulations(sim_setting, parameters = parameters_list[[i]], repetitions)
    sim_res$parameter <- param_iter[[i]]
    return(sim_res)
  })
  return(Reduce(rbind, df_res))
}



run_all_methods <- function(Adj_list, truecoms) {
  K <- length(unique(truecoms))
  
  # Note: to run graph-tool, uncomment the following lines and comment the next
  # methods_to_run <- c("dcmase", "ave_spherical", "sq-bias-adjusted",
  #                    "mase-spherical","lmfo", "graph-tool")
  # methods_to_run <- c("dcmase", "ave_spherical", "sq-bias-adjusted",
  #                     "mase-spherical","lmfo")
  methods_to_run <- c("frost-mf","frost-us","us","mf","lmfo","dcmase")
  results <- sapply(methods_to_run, function(method) {
    print(method)
    if (method=="frost-mf"){
      
      
      seed_for_python <- sample.int(.Machine$integer.max, 1)
      
      Adj_list_numpy <- lapply(Adj_list, function(X) np$array(as.matrix(X)))
      X_list <- r_to_py(Adj_list_numpy)
      truecoms_np <- np$array(truecoms)
      
      
      res <- frost$frost_multilayer(X_list, K, init_method='MF-SC-CA',init_seed=seed_for_python)
      labels <- py_to_r(res[[2]])
      labels <- labels + 1
      
      # Sauvegarder X et truecoms dans un fichier .npz
      #np$savez("data_for_python.npz", X = X_np, truecoms = truecoms_np,moments=labels)
      
      
      classError(labels, truecoms)$errorRate
      
    }
    else if (method=="frost-us"){
      seed_for_python <- sample.int(.Machine$integer.max, 1)
      
      Adj_list_numpy <- lapply(Adj_list, function(X) np$array(as.matrix(X)))
      X_list <- r_to_py(Adj_list_numpy)
      truecoms_np <- np$array(truecoms)
      
      
      res <- frost$frost_multilayer(X_list, K, init_method='USENC',init_seed=seed_for_python)
      labels <- py_to_r(res[[2]])
      labels <- labels + 1
      
      # Sauvegarder X et truecoms dans un fichier .npz
      #np$savez("data_for_python.npz", X = X_np, truecoms = truecoms_np,moments=labels)
      
      
      classError(labels, truecoms)$errorRate
      
    } 
    else if (method=="mf"){
     
      
      seed_for_python <- sample.int(.Machine$integer.max, 1)
      
      Adj_list_numpy <- lapply(Adj_list, function(X) np$array(as.matrix(X)))
      X_list <- r_to_py(Adj_list_numpy)
      truecoms_np <- np$array(truecoms)
      
      
      res <- frost$frost_multilayer(X_list, K, maxiter=as.integer(0), init_method='MF-SC-CA',init_seed=seed_for_python)
      labels <- py_to_r(res[[2]])
      labels <- labels + 1
      
      # Sauvegarder X et truecoms dans un fichier .npz
      #np$savez("data_for_python.npz", X = X_np, truecoms = truecoms_np,moments=labels)
      
      
      classError(labels, truecoms)$errorRate
      
      
    } 
    else if (method=="us"){
      
      
      seed_for_python <- sample.int(.Machine$integer.max, 1)
      
      Adj_list_numpy <- lapply(Adj_list, function(X) np$array(as.matrix(X)))
      X_list <- r_to_py(Adj_list_numpy)
      truecoms_np <- np$array(truecoms)
      
      
      res <- frost$frost_multilayer(X_list, K, maxiter=as.integer(0),init_method='USENC',init_seed=seed_for_python)
      labels <- py_to_r(res[[2]])
      labels <- labels + 1
      
      # Sauvegarder X et truecoms dans un fichier .npz
      #np$savez("data_for_python.npz", X = X_np, truecoms = truecoms_np,moments=labels)
      
      
      classError(labels, truecoms)$errorRate
      
    } 
    else{
     
      classError(comdetmethods(Adj_list, K, method = method), truecoms)$errorRate
     
    }
    
    
  }
  )
  # Note: to run graph-tool, uncomment the following lines and comment the next
  #names(results) <- c("DC-MASE", "Sum A", "S-A^2-Bias-adj",
  #                    "MASE", "OLMF", "graph-tool")
  
  names(results) <-  c("FROST_MF","FROST_US","US","MF","LMFO","DC_MASE")
  return(results)
}



#########################################################################
### Simulation 1: same B same theta
simulation1 <- function(parameters, seed = 1989) {
  set.seed(seed)
  m <- parameters[1]
  ave_deg <- 10
  
  # Generate network parameters ------------------------------------------------
  n <- 150
  K <- 3
  Z <- kronecker(diag(K), rep(1, n/K))
  #theta <- runif(n, min = 0.05, max = 1)
  theta <- rexp(n) + 0.2
  theta <- as.vector(theta / (Z%*%crossprod(Z,theta)/(n/K)))
  
  B <- 0.06*diag(K) + 0.04
  
  degree_corrections <- lapply(1:m, function(i) theta)
  B_matrices <- lapply(1:m, function(i) B)
  
  # Sample graphs --------------------------------------------------------------
  Adj_list <- lapply(1:m, function(i) {
    P <- tcrossprod((degree_corrections[[i]] * Z) %*% B_matrices[[i]], degree_corrections[[i]] * Z)
    P <- (ave_deg*n/sum(P)) * P
    sample_from_P(P)
  })
  return(list(Adj_list = Adj_list, truecom = as.vector(Z %*% 1:K)))
}

### Simulation 2: same theta, different B
simulation2 <- function(parameters, seed = 1989) {
  set.seed(seed)
  m <- parameters[1]
  ave_deg <- 10
  # Generate network parameters ------------------------------------------------
  n <- 150
  K <- 3
  Z <- kronecker(diag(K), rep(1, n/K))
  
  #theta <- runif(n, min = 0.05, max = 1)
  theta <- rexp(n) + 0.2
  theta <- as.vector(theta / (Z%*%crossprod(Z,theta)/(n/K)))
  degree_corrections <- lapply(1:m, function(i) theta)
  
  B_matrices <- lapply(1:m, function(i) {
    p <- runif(1)
    q <- runif(1)
    (p-q) *diag(K) + q
  })
  
  # Sample graphs --------------------------------------------------------------
  Adj_list <- lapply(1:m, function(i) {
    P <- tcrossprod((degree_corrections[[i]] * Z) %*% B_matrices[[i]], degree_corrections[[i]] * Z)
    P <- (ave_deg*n/sum(P)) * P
    sample_from_P(P)
  })
  return(list(Adj_list = Adj_list, truecom = as.vector(Z %*% 1:K)))
}





### Simulation 3: different B, different theta
simulation3 <- function(parameters, seed = 1989) {
  set.seed(seed)
  m <- parameters[1]
  ave_deg <- 10
  # Generate network parameters ------------------------------------------------
  n <- 150
  K <- 3
  Z <- kronecker(diag(K), rep(1, n/K))
  #degree_corrections <- lapply(1:m, function(i) runif(n, min = sqrt(0.05), max = 1)^2)
  degree_corrections <- lapply(1:m, function(i) {
    theta <- runif(n, min = 0.05, max = 1)
    theta <- rexp(n) + 0.2
    as.vector(theta / (Z%*%crossprod(Z,theta)/(n/K)))
  })
  B_matrices <- lapply(1:m, function(i) {
    p <- runif(1)
    q <- runif(1)
    (p-q) *diag(K) + q
  })
  
  # Sample graphs --------------------------------------------------------------
  Adj_list <- lapply(1:m, function(i) {
    P <- tcrossprod((degree_corrections[[i]] * Z) %*% B_matrices[[i]], degree_corrections[[i]] * Z)
    P <- (ave_deg*n/sum(P)) * P
    sample_from_P(P)
  })
  return(list(Adj_list = Adj_list, truecom = as.vector(Z %*% 1:K)))
}


### Simulation 4: same B, different theta
simulation4 <- function(parameters, seed = 1989) {
  set.seed(seed)
  m <- parameters[1]
  ave_deg <- 10
  # Generate network parameters ------------------------------------------------
  n <- 150
  K <- 3
  Z <- kronecker(diag(K), rep(1, n/K))
  
  degree_corrections <- lapply(1:m, function(i) {
    theta <- runif(n, min = 0.05, max = 1)
    theta <- rexp(n) + 0.2
    as.vector(theta / (Z%*%crossprod(Z,theta)/(n/K)))
  })
  B_matrices <- lapply(1:m, function(i) {
    B <- 0.06*diag(K) + 0.04
    B
  })
  
  # Sample graphs --------------------------------------------------------------
  Adj_list <- lapply(1:m, function(i) {
    P <- tcrossprod((degree_corrections[[i]] * Z) %*% B_matrices[[i]], degree_corrections[[i]] * Z)
    P <- (ave_deg*n/sum(P)) * P
    sample_from_P(P)
  })
  return(list(Adj_list = Adj_list, truecom = as.vector(Z %*% 1:K)))
}






# Simulation 6
### Simulation 6: changing high degree
simulation6 <- function(parameters, seed = 1989) {
  set.seed(seed)
  m <- parameters[1]
  # Set network parameters ------------------------------------------------
  n <- parameters[2] #150
  K <- parameters[3] #3
  
  ave_deg <- 10
  Z <- kronecker(diag(K), rep(1, n/K))
  
  theta1 <- rep(c(rep(0.15, 0.5*n/K), rep(0.8, 0.5*n/K)), K)
  theta2 <- 0.95 - theta1
  theta1 <- as.vector(theta1 / (Z%*%crossprod(Z,theta1)/(n/K)))
  theta2 <- as.vector(theta2 / (Z%*%crossprod(Z,theta2)/(n/K)))               
  degree_corrections <- lapply(1:m, function(i) if(i%%2==1){
    theta1
  }else{
    theta2
  })
  
  B_matrices <- lapply(1:m, function(i) {
    B <- 0.06*diag(K) + 0.04
    B
  })
  
  
  # Sample graphs --------------------------------------------------------------
  Adj_list <- lapply(1:m, function(i) {
    P <- tcrossprod((degree_corrections[[i]] * Z) %*% B_matrices[[i]], degree_corrections[[i]] * Z)
    P <- (ave_deg*n/sum(P)) * P
    sample_from_P(P)
  })
  #plot_adjmatrix(Adj_list[[1]])
  #colSums(Adj_list[[1]])
  truecom <- as.vector(Z %*% 1:K)
  
  results <- list(Adj_list = Adj_list, truecom = truecom)
  return(results)
}





# Simulation 7
### Simulation 7: changing high degree and random B
simulation7 <- function(parameters, seed = 1989) {
  set.seed(seed)
  m <- parameters[1]
  # Set network parameters ------------------------------------------------
  n <- parameters[2] #150
  K <- parameters[3] #3
  
  ave_deg <- 10
  Z <- kronecker(diag(K), rep(1, n/K))
  
  theta1 <- rep(c(rep(0.15, 0.5*n/K), rep(0.8, 0.5*n/K)), K)
  theta2 <- 0.95 - theta1
  theta1 <- as.vector(theta1 / (Z%*%crossprod(Z,theta1)/(n/K)))
  theta2 <- as.vector(theta2 / (Z%*%crossprod(Z,theta2)/(n/K)))               
  degree_corrections <- lapply(1:m, function(i) if(i%%2==1){
    theta1
  }else{
    theta2
  })
  
  B_matrices <- lapply(1:m, function(i) {
    p <- runif(1)
    q <- runif(1)
    ((p-q) *diag(K) + q)
  })
  
  # Sample graphs --------------------------------------------------------------
  Adj_list <- lapply(1:m, function(i) {
    P <- tcrossprod((degree_corrections[[i]] * Z) %*% B_matrices[[i]], degree_corrections[[i]] * Z)
    P <- (ave_deg*n/sum(P)) * P
    sample_from_P(P)
  })
  #plot_adjmatrix(Adj_list[[1]])
  #colSums(Adj_list[[1]])
  truecom <- as.vector(Z %*% 1:K)
  
  results <- list(Adj_list = Adj_list, truecom = truecom)
  return(results)
}