library(Rcpp)
library(MASS)
library(invgamma)
library(truncnorm)
library(fBasics)
library(Rcpp)
library(coda)
library(lme4)
library(binaryLogic)
library(BartRF)

lower_matrix <- function(a) {
  n <- length(a)
  dimen <- (1+sqrt(1+8*n))/2
  lower_triangle_matrix <- matrix(0, nrow = dimen, ncol = dimen)
  diag(lower_triangle_matrix) <- 1
  for (i in 2:dimen) {
    lower_triangle_matrix[i,1:(i-1)] <- a[(1+(i-1)*(i-2)/2):(1+(i-1)*(i-2)/2+(i-2))]
  }
  return(lower_triangle_matrix)
}

# Load necessary library
library(MASS) # for mvrnorm function
library(nnet)

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
  # correlation_factor <- 0.2
  # X_all[, 6] <- correlation_factor * Z_numeric[, 1] + rnorm(N, 0, sqrt(1 - correlation_factor^2))  # X6 ~ Z1
  

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
  
  # Z_useful <- cbind(as.matrix(Z_numeric[, 1:useful_z]), 1)
  # useful_number <- ncol(Z_useful)
  # beta_matrix <- matrix(0, nrow = N, ncol = useful_number)
  # 
  # beta_correlation_factor <- 0.4  # smaller to limit beta-Z association
  # 
  # for (j in 1:useful_number) {
  #   for (i in 1:n) {
  #     if (j == useful_number) {
  #       cluster_beta <- rnorm(1, mean = 0, sd = 0.1)
  #     } else {
  #       # Z_val <- Z_cluster[i, j]
  #       # shrink_factor <- 0.4  # <-- tune this to shrink variance, preserves correlation
  #       # Z_val_shrunk <- Z_val * shrink_factor
  #       # cluster_beta <- rnorm(1, mean = beta_correlation_factor * Z_val_shrunk, sd = 0.2)
  # 
  #       Z_val <- Z_cluster[i, j]
  #       cluster_beta <- rnorm(1, mean = beta_correlation_factor * Z_val, sd = 0.2)
  # 
  # 
  #     }
  #     beta_matrix[(m * (i - 1) + 1):(m * i), j] <- cluster_beta
  #   }
  # }

  
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
  y <- scale(fixed_effect + random_effect / 3) + noise
  
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

taskid <- as.numeric(commandArgs(trailingOnly = TRUE))

print(taskid)

data <- generate_rf_test(taskid,50,100,3,4)

q = ncol(rfModelMatrix(data$Z)$Z)
c0 <- 0.05   # shape
d0 <- 0.05   # rate
## b: standard normal
## gamma: multivariate normal distribution
gamma_mean  <- rep(0, q*(q-1)/2)
gamma_cov <- 0.5 *diag(q*(q-1)/2)
## lambda: truncated normal distribution
z_alpha <- rep(0.5,q)
z_beta <- rep(0.5,q)
lambda_mean <- rep(0, q)
lambda_cov <- rep(30, q) #should be square

test_Z <- data$Z
test_Z[,1:2] <- scale(test_Z[,1:2])
test_Z[,4:7] <- scale(test_Z[,4:7])

result <- BartRF::permute.vs(x.train = as.data.frame(data$X), y.train = data$y, z.train = as.data.frame(test_Z), id = data$bart_id, nreps = 10, npermute = 100,
                             nskip = 5000, ndpost = 5000, ntree = 20,
                             z_c0 = c0, z_d0 = d0,
                             z_gamma_mean = gamma_mean, z_gamma_cov = gamma_cov,
                             z_lambda_mean = lambda_mean, z_lambda_cov = lambda_cov, z_alpha = z_alpha, z_beta = z_beta, plot = FALSE)

#save result
saveRDS(result, file = paste0("10_4/result_",taskid,".rds"))
