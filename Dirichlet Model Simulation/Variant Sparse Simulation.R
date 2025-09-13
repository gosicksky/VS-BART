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
library(LaplacesDemon)

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
  #sigma <- rinvwishart(useful_number + 2, diag(useful_number))
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
  y <- fixed_effect + random_effect / 3 + noise
  
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
taskid <- as.numeric(commandArgs(trailingOnly = TRUE))

print(taskid)

data <- generate_rf_test(taskid,50,100,3,3)

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
test_Z[,4:6] <- scale(test_Z[,4:6])

result <- BartRF::wbart(x.train = as.data.frame(data$X), y.train = data$y, z.train = as.data.frame(test_Z), id = data$bart_id,
                             nskip = 5000, ndpost = 5000, ntree = 50,
                             z_c0 = c0, z_d0 = d0,
                             z_gamma_mean = gamma_mean, z_gamma_cov = gamma_cov,
                             z_lambda_mean = lambda_mean, z_lambda_cov = lambda_cov, z_alpha = z_alpha, z_beta = z_beta, sparse = TRUE)
result$Sigma <- data$Sigma

result1 <- list(varprob = result$varprob, lambda = result$combined.lambda)

#save result
saveRDS(result1, file = paste0("10_3_sparse_50/result_",taskid,".rds"))
