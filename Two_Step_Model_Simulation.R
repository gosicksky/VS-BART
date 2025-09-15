library(Rcpp)
library(MASS)
library(invgamma)
library(truncnorm)
library(fBasics)
library(coda)
library(lme4)
library(binaryLogic)
library(Matrix)
library(BH)
library(RcppEigen)
library(BartMixVs)

sourceCpp("./src/cwbart.cpp", rebuild = FALSE, verbose = TRUE)
source("./R/wbart.R")
source("./R/rfModleMatrix.R")
source("./R/bartModelMatrix.R")
source("./R/permute.vs.R")

library(nnet)
library(glmnet)

dic = "./simulation_results/"
args = commandArgs(trailingOnly = TRUE)
iter = as.integer(args[1])

generate_rf_test <- function(seed, cluster_number, cluster_size, useful_z, noise_z) {
  set.seed(seed)
  
  n <- cluster_number
  m <- cluster_size
  N <- n * m
  rf_number <- useful_z + noise_z
  p <- 10
  bart_id <- rep(1:n, each = m)
  
  f0 <- function(x) {
    term1 <- - ((as.numeric(x[,7])) ^ 2) + 3 * (as.numeric(x[,5]))
    term2 <- 6 * (as.numeric(x[,1]) - 1)
    term3 <- as.numeric(x[,6]) * as.numeric(x[,8])
    return(term1 + term2 + term3)
  }
  
  
  useful_x_indices <- c(1, 5, 6, 7, 8)
  noise_x_indices <- setdiff(1:p, useful_x_indices)
  
  # Step 1: Generate cluster-level Z
  Z_cluster <- matrix(0, nrow = n, ncol = rf_number)
  for (i in 1:useful_z) {
    if (i == useful_z) {
      Z_cluster[, i] <- rbinom(n, 3, 0.5)
    } else {
      Z_cluster[, i] <- rnorm(n, 0, 1)
    }
  }
  
  for (i in (useful_z + 1):rf_number) {
    Z_cluster[, i] <- rnorm(n, 0, 1)
  }
  
  # Step 2: Expand Z to individual level
  Z_numeric <- Z_cluster[bart_id, ]
  Z <- as.data.frame(Z_numeric)
  Z[, useful_z] <- as.factor(Z[, useful_z])
  
  # Step 3: Generate individual-level X
  
  X_all <- matrix(rnorm(N * p, 0, 1), nrow = N, ncol = p)
  
  correlation_factor <- 0.2
  X_all[, 5] <- correlation_factor * Z_numeric[, 1] + rnorm(N, 0, sqrt(1 - correlation_factor^2))  # X5 ~ Z1
  X_all[, 6] <- correlation_factor * Z_numeric[, 2] + rnorm(N, 0, sqrt(1 - correlation_factor^2))  # X6 ~ Z2
  
  X <- as.data.frame(X_all)
  binary_positions <- useful_x_indices[1]
  for (i in binary_positions) {
    X[, i] <- as.factor(rbinom(N, 1, 0.5))
  }
  
  # Step 4: Generate beta_matrix using small variance and weak correlation with Z
  library(MASS)
  
  Z_useful <- cbind(as.matrix(Z_numeric[, 1:useful_z]), 1)  # add intercept column
  useful_number <- ncol(Z_useful)  # includes intercept
  beta_matrix <- matrix(0, nrow = N, ncol = useful_number)
  
  rho_beta <- -0.8
  K <- useful_z
  
  R <- matrix(rho_beta, nrow = K, ncol = K)
  diag(R) <- 1
  is_pos_def <- function(mat) {
    all(eigen(mat, symmetric = TRUE, only.values = TRUE)$values > 0)
  }
  
  while (!is_pos_def(R)) {
    rho_beta <- rho_beta + 0.05
    R <- matrix(rho_beta, nrow = K, ncol = K)
    diag(R) <- 1
  }
  custom_vars <- c(0.01, 0.01, 0.01) 
  scaling_matrix <- diag(sqrt(custom_vars))
  cov_beta <- scaling_matrix %*% R %*% scaling_matrix  # Σ = D R D
  
  beta_by_cluster <- matrix(NA, nrow = n, ncol = useful_number)
  
  beta_correlation_factors <- c(2.5, 2.5, 2)
  
  for (i in 1:n) {
    Z_vals <- Z_cluster[i, 1:useful_z]
    beta_mean_slopes <- beta_correlation_factors * Z_vals
    beta_slopes <- mvrnorm(1, mu = beta_mean_slopes, Sigma = cov_beta)
    beta_intercept <- rnorm(1, mean = 0, sd = 0.1)
    beta_by_cluster[i, ] <- c(beta_slopes, beta_intercept)
  }
  
  beta_matrix <- beta_by_cluster[bart_id, ]
  
  
  
  # Step 5: Final outcome
  noise <- rnorm(N, 0, 1)
  fixed_effect <- f0(X)
  random_effect <- rowSums(beta_matrix * Z_useful)
  y <- scale(fixed_effect + random_effect + noise)
  
  rf_true_variable <- 1:useful_z
  Z[, useful_z] <- as.factor(Z[, useful_z])
  
  return(list(
    X = X,
    y = y,
    bart_id = bart_id,
    Z = Z,
    beta_matrix = beta_matrix,
    noise = noise,
    random_effect = random_effect,
    fix_effect = fixed_effect,
    rf_true_variable = rf_true_variable
  ))
}

# Define the sensitivity analysis alpha values
alpha_values <- c(0.005, 0.01, 0.05, 0.1)
ntree_values <- c(20)

# Function to run a simulation with a specific alpha threshold
run_simulation <- function(
    dic,
    iter,
    cluster_number, 
    cluster_size, 
    useful_z, 
    noise_z, 
    ntree,
    alpha,  # Added parameter for sensitivity analysis
    c0 = 0.05, 
    d0 = 0.05, 
    z_alpha = rep(0.5, 1), 
    z_beta = rep(0.5, 1), 
    lambda_mean = rep(0, 1), 
    lambda_cov = rep(30, 1), 
    gamma_mean = rep(0, 1), 
    gamma_cov = 1*diag(1)
) {
  
  if (!dir.exists(dic)) {
    dir.create(dic, recursive = TRUE)
  }
  
  library(BartMixVs)
  library(BartRI)
  
  total_variables <- useful_z + noise_z
  vip_counts_indep <- numeric(total_variables)
  mi_counts_indep <- numeric(total_variables)
  within_vip_counts_indep <- numeric(total_variables)  
  
  set.seed(iter)
  data_rf_test <- generate_rf_test(
    seed = iter,
    cluster_number = cluster_number,
    cluster_size = cluster_size,
    useful_z = useful_z,
    noise_z = noise_z
  )
  
  set.seed(iter)
  test_Z_indep <- matrix(rep(1, nrow(data_rf_test$Z)))
  
  bart_rf_test_indep <- BartRI::wbart(
    x.train = as.data.frame(data_rf_test$X), 
    y.train = data_rf_test$y, 
    nskip = 15000, 
    ndpost = 5000, 
    ntree = ntree,
    z.train = cbind(test_Z_indep), 
    id = data_rf_test$bart_id,
    z_c0 = c0, 
    z_d0 = d0,
    z_gamma_mean = gamma_mean, 
    z_gamma_cov = gamma_cov,
    z_lambda_mean = lambda_mean, 
    z_lambda_cov = lambda_cov, 
    z_alpha = z_alpha, 
    z_beta = z_beta, 
    rm.const = FALSE
  )
  
  b_indep <- colMeans(bart_rf_test_indep$rf_b * matrix(bart_rf_test_indep$original.lambda, nrow = 5000, ncol = cluster_number))
  z.train_indep <- as.data.frame(unique(data_rf_test$Z))
  
  # Apply variable selection with the current alpha threshold
  permute_indep <- BartMixVs::permute.vs(
    x.train = z.train_indep, 
    y.train = b_indep, 
    ndpost = 5000, 
    nskip = 5000,
    alpha = alpha  # Using variable alpha for sensitivity analysis
  )
  
  vip_selected_indep <- as.integer(gsub("V", "", permute_indep$vip.imp.names))  
  mi_selected_indep <- as.integer(gsub("V", "", permute_indep$mi.imp.names))    
  within_selected_indep <- as.integer(gsub("V", "", permute_indep$within.type.vip.imp.names))
  
  # Initialize counts for each method
  vip_counts_indep[vip_selected_indep] <- vip_counts_indep[vip_selected_indep] + 1
  mi_counts_indep[mi_selected_indep] <- mi_counts_indep[mi_selected_indep] + 1
  within_vip_counts_indep[within_selected_indep] <- within_vip_counts_indep[within_selected_indep] + 1  # New for within-type VIP
  
  # Convert selected variables to strings
  selected_submodel_vip <- paste0("Z", vip_selected_indep, collapse = ",")
  selected_submodel_mi <- paste0("Z", mi_selected_indep, collapse = ",")
  selected_submodel_within <- paste0("Z", within_selected_indep, collapse = ",")
  
  # If want to output the X selection result
  # permute_indep_x <- BartRI::permute.vs(
  #     x.train = as.data.frame(data_rf_test$X), 
  #     y.train = data_rf_test$y,
  #     z.train = cbind(test_Z_indep), 
  #     ndpost = 5000, 
  #     nskip = 5000,
  #     alpha = 0.05,
  #    id = data_rf_test$bart_id,
  #     z_c0 = c0, 
  #     z_d0 = d0,
  #     z_gamma_mean = gamma_mean, 
  #     z_gamma_cov = gamma_cov,
  #     z_lambda_mean = lambda_mean, 
  #     z_lambda_cov = lambda_cov, 
  #     z_alpha = z_alpha, 
  #     z_beta = z_beta
  #   )
  # 
  #   x_vip_selected <- as.integer(gsub("V", "", permute_indep_x$vip.imp.names))
  #   x_mi_selected <- as.integer(gsub("V", "", permute_indep_x$mi.imp.names))
  #   x_within_selected <- as.integer(gsub("V", "", permute_indep_x$within.type.vip.imp.names))
  # 
  # x_vip_counts <- numeric(ncol(data_rf_test$X))
  # x_mi_counts <- numeric(ncol(data_rf_test$X))
  # x_within_vip_counts <- numeric(ncol(data_rf_test$X))
  # 
  # x_vip_counts[x_vip_selected] <- x_vip_counts[x_vip_selected] + 1
  # x_mi_counts[x_mi_selected] <- x_mi_counts[x_mi_selected] + 1
  # x_within_vip_counts[x_within_selected] <- x_within_vip_counts[x_within_selected] + 1
  # 
  #   selected_submodel_x_vip <- paste0("X", x_vip_selected, collapse = ",")
  #   selected_submodel_x_mi <- paste0("X", x_mi_selected, collapse = ",")
  #   selected_submodel_x_within <- paste0("X", x_within_selected, collapse = ",")
  
  # Define function to compute evaluation metrics
  compute_metrics <- function(selected_vars, true_vars) {
    tp <- length(intersect(selected_vars, true_vars))  # True Positives
    fp <- length(setdiff(selected_vars, true_vars))    # False Positives
    fn <- length(setdiff(true_vars, selected_vars))    # False Negatives
    
    precision <- ifelse(tp + fp > 0, tp / (tp + fp), 0)
    recall <- ifelse(tp + fn > 0, tp / (tp + fn), 0)
    f1_score <- ifelse(precision + recall > 0, 2 * precision * recall / (precision + recall), 0)
    type_II_error <- ifelse(tp + fn > 0, fn / (tp + fn), 0)
    error_rate <- ifelse(tp + fp + fn > 0, (fp + fn) / (tp + fp + fn), 0)
    type_I_error <- ifelse(noise_z > 0, fp / noise_z, NA)  
    
    return(list(type_I = type_I_error, precision = precision, recall = recall, f1 = f1_score, type_II = type_II_error, error = error_rate))
  }
  
  # Compute metrics for each method
  metrics_vip <- compute_metrics(vip_selected_indep, data_rf_test$rf_true_variable)
  metrics_mi <- compute_metrics(mi_selected_indep, data_rf_test$rf_true_variable)
  metrics_within <- compute_metrics(within_selected_indep, data_rf_test$rf_true_variable)
  
  # true_fixed_vars <- c(1, 5, 6, 7, 8)
  #   metrics_x_vip <- compute_metrics(x_vip_selected, true_fixed_vars)
  #   metrics_x_mi <- compute_metrics(x_mi_selected, true_fixed_vars)
  #   metrics_x_within <- compute_metrics(x_within_selected, true_fixed_vars)
  
  # Write results to CSV
  write.csv(x = data.frame(
    ntree = ntree,
    alpha = alpha,
    vip_counts_indep = vip_counts_indep, 
    mi_counts_indep = mi_counts_indep, 
    within_vip_counts_indep = within_vip_counts_indep,
    
    selected_submodel_vip = selected_submodel_vip,
    selected_submodel_mi = selected_submodel_mi, 
    selected_submodel_within = selected_submodel_within,
    
    precision_vip = metrics_vip$precision,
    recall_vip = metrics_vip$recall,
    f1_vip = metrics_vip$f1,
    type_II_vip = metrics_vip$type_II,
    error_rate_vip = metrics_vip$error,
    type_I_vip = metrics_vip$type_I,
    
    precision_mi = metrics_mi$precision,
    recall_mi = metrics_mi$recall,
    f1_mi = metrics_mi$f1,
    type_II_mi = metrics_mi$type_II,
    error_rate_mi = metrics_mi$error,
    type_I_mi = metrics_mi$type_I,
    
    precision_within = metrics_within$precision,
    recall_within = metrics_within$recall,
    f1_within = metrics_within$f1,
    type_II_within = metrics_within$type_II,
    error_rate_within = metrics_within$error,
    type_I_within = metrics_within$type_I
    
    
  ), paste0(dic, "/iter_", iter, "_alpha_", alpha,"_ntree_", ntree, "_cluster_number_", cluster_number, "_cluster_size_", cluster_size, 
            "_useful_z_", useful_z, "_noise_z_", noise_z, ".csv"))
  
  
}


# x_vip_counts = paste(x_vip_counts, collapse = ","),
# x_mi_counts = paste(x_mi_counts, collapse = ","),
# x_within_vip_counts = paste(x_within_vip_counts, collapse = ","),
# 
# selected_submodel_x_vip = selected_submodel_x_vip,
# selected_submodel_x_mi = selected_submodel_x_mi,
# selected_submodel_x_within = selected_submodel_x_within,
#     
# precision_x_vip = metrics_x_vip$precision,
# recall_x_vip = metrics_x_vip$recall,
# f1_x_vip = metrics_x_vip$f1,
# type_II_x_vip = metrics_x_vip$type_II,
# error_rate_x_vip = metrics_x_vip$error,
#     
# precision_x_mi = metrics_x_mi$precision,
# recall_x_mi = metrics_x_mi$recall,
# f1_x_mi = metrics_x_mi$f1,
# type_II_x_mi = metrics_x_mi$type_II,
# error_rate_x_mi = metrics_x_mi$error,
#     
# precision_x_within = metrics_x_within$precision,
# recall_x_within = metrics_x_within$recall,
# f1_x_within = metrics_x_within$f1,
# type_II_x_within = metrics_x_within$type_II,
# error_rate_x_within = metrics_x_within$error


# Loop over different alpha values for sensitivity analysis
for (ntree in ntree_values) {
  for (alpha in alpha_values) {
    run_simulation(
      dic = paste0(dic, "/sim1/"),
      iter = iter,
      cluster_number = 50, 
      cluster_size = 100, 
      useful_z = 3, 
      noise_z = 3, 
      alpha = alpha,
      ntree = ntree
    )
    
    run_simulation(
      dic = paste0(dic, "/sim2/"),
      iter = iter,
      cluster_number = 50, 
      cluster_size = 100, 
      useful_z = 3, 
      noise_z = 4, 
      alpha = alpha,
      ntree = ntree
    )
    
    run_simulation(
      dic = paste0(dic, "/sim3/"),
      iter = iter,
      cluster_number = 50, 
      cluster_size = 100, 
      useful_z = 3, 
      noise_z = 5, 
      alpha = alpha,
      ntree = ntree
    )
    
    run_simulation(
      dic = paste0(dic, "/sim4/"),
      iter = iter,
      cluster_number = 50, 
      cluster_size = 100, 
      useful_z = 3, 
      noise_z = 6, 
      alpha = alpha,
      ntree = ntree
    )
  }
}