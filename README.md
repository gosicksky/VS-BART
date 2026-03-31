Multi-level variable selection using a BART-enhanced mixed-effects framework
================================================================================

This is the repository for the paper "Multi-level variable selection using a BART-enhanced mixed-effects framework".

## Install the package

```
install.packages("BartVS_0.1.0.tar.gz")
```

BartRF and BartRI is only for the code submitted in revision.

## Functions

The 'BartVS' package provides the following functions for variable selection in multi-level models:

1. `permute.vs_unified` implements the variable selection procedure for the unified model, which simultaneously models fixed and random effects.

2. `permute.vs_two_step` implements the variable selection procedure for the two-step model, which first estimates random effects and then performs variable selection on the fitted values.

3. `wabrt_unified` is a wrapper function for fitting the unified BART model, which is used within the variable selection procedure. The parameter `sparse` controls whether to use a sparse BART model for variable selection.

4. `wbart_two_step` is a wrapper function for fitting the two-step BART model, which is used within the variable selection procedure for the two-step model and could also be used for random effect variable selection as the second step.

## Example

Data example can be seen in data_simulation.R. The following code snippets show how to use the functions for variable selection in the unified model, two-step model, and sparse BART variable selection.

##### Unified model variable selection

```R
q = ncol(BartVS::rfModelMatrix_unified(data$Z)$Z)
## gamma: multivariate normal distribution
gamma_mean  <- rep(0, q*(q-1)/2)
gamma_cov <- 0.5 *diag(q*(q-1)/2)
## lambda: truncated normal distribution
z_alpha <- rep(0.5,q)
z_beta <- rep(0.5,q)
lambda_mean <- rep(0, q)
lambda_cov <- rep(30, q)
result <- BartVS::permute.vs_unified(x.train = as.data.frame(data$X), y.train = data$y, z.train = as.data.frame(data$Z), id = data$id, nreps = 10, npermute = 100,
                             nskip = 5000, ndpost = 5000, ntree = 20,
                             z_gamma_mean = gamma_mean, z_gamma_cov = gamma_cov,
                             z_lambda_mean = lambda_mean, z_lambda_cov = lambda_cov, z_alpha = z_alpha, z_beta = z_beta)

```

##### Two-step model variable selection

```R
z_intercept <- matrix(rep(1, nrow(data$Z)))
## lambda: truncated normal distribution
z_alpha = rep(0.5, 1)
z_beta = rep(0.5, 1)
lambda_mean = rep(0, 1)
lambda_cov = rep(30, 1)
#permutation test for variable selection
fix_result <- BartVS::permute.vs_two_step(x.train = as.data.frame(data$X), y.train = data$y, z.train = as.data.frame(z_intercept), id = data$id, nreps = 10, npermute = 100,
                             nskip = 5000, ndpost = 5000, ntree = 20,
                             z_lambda_mean = lambda_mean, z_lambda_cov = lambda_cov, z_alpha = z_alpha, z_beta = z_beta)
#second step for random effects variable selection
random_result <- BartVS::wbart_two_step(x.train = as.data.frame(data$X), y.train = data$y, z.train = as.data.frame(z_intercept), id = data$id, nskip = 5000, ndpost = 5000, ntree = 20,
                             z_lambda_mean = lambda_mean, z_lambda_cov = lambda_cov, z_alpha = z_alpha, z_beta = z_beta, two_step_alpha = 0.05)
```

##### Sparse BART variable selection

```R
q = ncol(BartVS::rfModelMatrix_unified(data$Z)$Z)
## gamma: multivariate normal distribution
gamma_mean  <- rep(0, q*(q-1)/2)
gamma_cov <- 0.5 *diag(q*(q-1)/2)
## lambda: truncated normal distribution
z_alpha <- rep(0.5,q)
z_beta <- rep(0.5,q)
lambda_mean <- rep(0, q)
lambda_cov <- rep(30, q)
result <- BartVS::wbart_unified(x.train = as.data.frame(data$X), y.train = data$y, z.train = as.data.frame(data$Z), id = data$id, nskip = 5000, ndpost = 5000, ntree = 20,
                             z_gamma_mean = gamma_mean, z_gamma_cov = gamma_cov,
                             z_lambda_mean = lambda_mean, z_lambda_cov = lambda_cov, z_alpha = z_alpha, z_beta = z_beta,
                             sparse = TRUE)
```

## Future plans

Simplify the code and Upload the package to CRAN.