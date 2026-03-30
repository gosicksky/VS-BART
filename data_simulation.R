generate_unified <- function(seed, cluster_number, cluster_size, useful_z, noise_z) {
  set.seed(seed)
  
  n <- cluster_number
  m <- cluster_size
  N <- n * m
  rf_number <- useful_z + noise_z
  p <- 10
  bart_id <- rep(1:n, each = m)
  
  f0 <- function(x) {
    term1 <- - as.numeric(x[,7]) + 3 * (as.numeric(x[,5])) + ((as.numeric(x[,7])) ^ 2)
    term2 <- 6 * (as.numeric(x[,1]) - 1)
    term3 <- as.numeric(x[,6]) + as.numeric(x[,6]) * as.numeric(x[,8])
    return(term1 + term2 + term3)
  }
  
  X_all <- matrix(rnorm(N * p, 0, 1), nrow = N, ncol = p)
  X <- as.data.frame(X_all)
  for (i in c(1,2)) {
    X[, i] <- as.factor(rbinom(N, 1, 0.5))
  }
  
  # Step 1: two columns is subset of Z and others are independent
  Z <- matrix(0, nrow = N, ncol = rf_number)
  for (i in 1:useful_z) {
    if (i == useful_z) {
      Z[, i] <- rbinom(N, 3, 0.5)
    } else {
      Z[, i] <- rnorm(N, 0, 1)
    }
  }
  
  for (i in (useful_z + 1):rf_number) {
    Z[, i] <- rnorm(N, 0, 1)
  }
  
  # replace two columns with X
  rho <- 0.2
  epsilon <- rnorm(N, 0, 1)
  Z[, 1] <- rho * scale(X[,5]) + sqrt(1 - rho^2) * epsilon
  epsilon <- rnorm(N, 0, 1)
  Z[, 2] <- rho * scale(X[,6]) + sqrt(1 - rho^2) * epsilon
  
  Z_useful <- Z[, 1:useful_z]
  Z_useful <- cbind(1, Z_useful)
  
  useful_number = ncol(Z_useful)
  mu <- rep(0, useful_number)
  sigma <- rinvwishart(useful_number + 2, diag(useful_number))
  repeat {
    sigma <- rinvwishart(useful_number + 2, diag(useful_number))
    if (all(diag(sigma) > 0.2)) break
  }
  beta <- mvrnorm(n, mu, sigma)
  beta_matrix <- matrix(0, nrow = N, ncol = useful_number)
  for (i in 1:(N/m)) {
    beta_matrix[(m * (i - 1) + 1):(m * i),] <- rep(beta[i,], each = m)
  }
  
  # Step 5: Final outcome
  noise <- rnorm(N, 0, 1)
  fixed_effect <- f0(X)
  random_effect <- rowSums(beta_matrix * Z_useful)
  y <- fixed_effect + random_effect + noise
  
  rf_true_variable <- 1:useful_z
  Z <- as.data.frame(Z)
  Z[, useful_z] <- as.factor(Z[, useful_z])
  
  return(list(
    X = X,
    y = y,
    bart_id = bart_id,
    Z = Z,
    beta_matrix = beta_matrix,
    noise = noise,
    Sigma = sigma,
    random_effect = random_effect,
    fix_effect = fixed_effect,
    rf_true_variable = rf_true_variable
  ))
}


generate_two_step <- function(seed, cluster_number, cluster_size, useful_z, noise_z) {
  
  n <- cluster_number
  m <- cluster_size
  N <- n * m
  rf_number <- useful_z + noise_z
  p <- 10
  bart_id <- rep(1:n, each = m)
  
  
  set.seed(seed)
  f0 <- function(x) {
    term1 <- - as.numeric(x[,7]) + 3 * (as.numeric(x[,5])) + ((as.numeric(x[,7])) ^ 2)
    term2 <- 6 * (as.numeric(x[,1]) - 1)
    term3 <- as.numeric(x[,6]) + as.numeric(x[,6]) * as.numeric(x[,8])  
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
  binary_positions <- c(1,2)
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
  custom_vars <- c(0.01, 0.01, 0.01)  # adjust if useful_z ≠ 3
  scaling_matrix <- diag(sqrt(custom_vars))
  cov_beta <- scaling_matrix %*% R %*% scaling_matrix
  
  beta_by_cluster <- matrix(NA, nrow = n, ncol = useful_number)
  
  beta_correlation_factors <- c(2.5, 2.5, 2)
  
  for (i in 1:n) {
    Z_vals <- Z_cluster[i, 1:useful_z]
    beta_mean_slopes <- beta_correlation_factors * Z_vals
    beta_slopes <- mvrnorm(1, mu = beta_mean_slopes, Sigma = cov_beta)
    beta_intercept <- rnorm(1, mean = 0, sd = 1)
    beta_by_cluster[i, ] <- c(beta_slopes, beta_intercept)
  }
  
  beta_matrix <- beta_by_cluster[bart_id, ]
  
  # Step 5: Final outcome
  noise <- rnorm(N, 0, 1)
  fixed_effect <- f0(X)
  random_effect <- rowSums(beta_matrix * Z_useful)
  y <- fixed_effect + random_effect + noise
  
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