Multi-level variable selection using a BART-enhanced mixed-effects framework
================================================================================

This is the repository for the paper "Multi-level variable selection using a BART-enhanced mixed-effects framework".

## Install the package

```
install.packages("BartVS_0.1.0.tar.gz")
```

## Functions

The 'BartVS' package provides the following functions for variable selection in multi-level models:

`permute.vs_unified` implements the variable selection procedure for the unified model, which simultaneously models fixed and random effects.
`permute.vs_two_step` implements the variable selection procedure for the two-step model, which first estimates random effects and then performs variable selection on the fitted values.
`wabrt_unified` is a wrapper function for fitting the unified BART model, which is used within the variable selection procedure. The parameter `sparse` controls whether to use a sparse BART model for variable selection.

## Example

##### Unified model variable selection

```R
result <- BartVS::permute.vs_unified(x.train = as.data.frame(data$X), y.train = data$y, z.train = as.data.frame(data$Z), id = data$id, nreps = 10, npermute = 100,
                             nskip = 5000, ndpost = 5000, ntree = 20,
                             z_c0 = c0, z_d0 = d0,
                             z_gamma_mean = gamma_mean, z_gamma_cov = gamma_cov,
                             z_lambda_mean = lambda_mean, z_lambda_cov = lambda_cov, z_alpha = z_alpha, z_beta = z_beta)

```

##### Two-step model variable selection

```R
result <- BartVS::permute.vs_two_step(x.train = as.data.frame(data$X), y.train = data$y, z.train = as.data.frame(data$Z), id = data$id, nreps = 10, npermute = 100,
                             nskip = 5000, ndpost = 5000, ntree = 20,
                             z_c0 = c0, z_d0 = d0,
                             z_gamma_mean = gamma_mean, z_gamma_cov = gamma_cov,
                             z_lambda_mean = lambda_mean, z_lambda_cov = lambda_cov, z_alpha = z_alpha, z_beta = z_beta)
```

##### Sparse BART variable selection

```R
result <- BartVS::wbart_unified(x.train = as.data.frame(data$X), y.train = data$y, z.train = as.data.frame(data$Z), id = data$id, nskip = 5000, ndpost = 5000, ntree = 20,
                             z_c0 = c0, z_d0 = d0,
                             z_gamma_mean = gamma_mean, z_gamma_cov = gamma_cov,
                             z_lambda_mean = lambda_mean, z_lambda_cov = lambda_cov, z_alpha = z_alpha, z_beta = z_beta,
                             sparse = TRUE)
```

