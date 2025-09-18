Multi-level variable selection using a BART-enhanced mixed-effects framework
================================================================================

This is the repository for the paper "Multi-level variable selection using a BART-enhanced mixed-effects framework".

To use the code, you could use the BartRF_0.1.tar.gz (the unified model) and BartRI_0.1.tar.gz (the two-step model) zip files to install the packages in R.

## Example

For the unified model:

```r
result <- BartRF::permute.vs(x.train = as.data.frame(data$X), y.train = data$y, z.train = as.data.frame(data$Z), id = data$id, nreps = 10, npermute = 100,
                             nskip = 5000, ndpost = 5000, ntree = 20,
                             z_c0 = c0, z_d0 = d0,
                             z_gamma_mean = gamma_mean, z_gamma_cov = gamma_cov,
                             z_lambda_mean = lambda_mean, z_lambda_cov = lambda_cov, z_alpha = z_alpha, z_beta = z_beta)
```

For the two-step model: 

```r
result <- BartRI::wbart(x.train = as.data.frame(data$X), y.train = data$y, z.train = as.data.frame(random_intercept), id = data$id, nreps = 10, npermute = 100,
                             nskip = 5000, ndpost = 5000, ntree = 20,
                             z_c0 = c0, z_d0 = d0,
                             z_gamma_mean = gamma_mean, z_gamma_cov = gamma_cov,
                             z_lambda_mean = lambda_mean, z_lambda_cov = lambda_cov, z_alpha = z_alpha, z_beta = z_beta)
                             
b_indep <- colMeans(result$rf_b * matrix(result$original.lambda, nrow = 5000, ncol = cluster_number))
z.train_indep <- as.data.frame(unique(data$Z))
  
# Apply variable selection with the current alpha threshold
permute_indep <- BartMixVs::permute.vs(
    x.train = z.train_indep, 
    y.train = b_indep, 
    ndpost = 5000, 
    nskip = 5000,
    alpha = alpha  # Using variable alpha for sensitivity analysis
  )
```