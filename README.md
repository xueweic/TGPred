
<!-- README.md is generated from README.Rmd. Please edit that file -->

# APGD v.0.1.0

<!-- badges: start -->

The Accelerated Proximal Gradient Descent (APGD) algorithm is a R package to solve the penalized regression models, including 

- **HuberNet**: Huber loss function along with Network-based penalty function;
- **HuberLasso**: Huber loss function along with Lasso penalty function;
- **HuberENET**: Huber loss function along with Elastic Net penalty function;
- **ENET**: Mean square error loss function along with Elastic Net penalty function;
- **Lasso**: Mean square error loss function along with Lasso penalty function;
- **Net**: Mean square error loss function along with Network-based penalty function.

## Installation

You can install the released version of APGD from Github with:

``` r
devtools::install_github("xueweic/APGD")
```

## Reference
Xuewei Cao<sup>+</sup>, Ling Zhang<sup>+</sup>, Kui Zhang, Sanzhen Liu, Qiuying Sha*, Hairong Wei*. HuberNet function for interfering target genes of regulatory genes using high-throughput gene expression data.

<sub><sup> <sup>+</sup> These authors have contributed equally to this work </sup></sub>

**Any questions**? xueweic_AT_mtu_DOT_edu, lingzhan_AT_mtu_DOT_edu

## Examples

### 1. Simulated data

**Step 1**: Construct the network structure from either Hierarchical Network or Barabasi-Albert Network in simulation studies.

- In Hierarchical Network, the number of genes must be the integer times 100.
- In Barabasi-Albert Network, the number of genes must be the integer times 10.

``` r
library(APGD)
N_genes = 200
Adj = ConstructNetwork(N_genes, "HN")
Adj = ConstructNetwork(N_genes, "BAN")
```

**Step 2**: Calculate Laplacian matrix and symmetric normalized Laplacian matrix from an adjacency matrix.

``` r
Sigma1 = GraphicalModel(Adj)
```

**Step 3**: Simulate y and X from a given network structure (Adjacency matrix and Laplacian matrix).


``` r
N_samples <- 300
res = SimulationData(N_sample,N_genes,Adj,Sigma1,"BAN", beta0 = 1)
y = res$y
X = res$X
```
or 

``` r
res = SimulationData(N_sample,N_genes,Adj,Sigma1,"HN", beta0 = 1)
```

### 2. Estimate Regression Coefficients

Calculate the estimated regression coefficients $\beta$(beta_hat) using one of methods solving by APGD or CVX for a given set of 

- **HuberNet**: Huber loss function along with Network-based penalty function.

``` r
lambda0 = 200
alpha0 = 0.5
beta_hat_APGD <- HuberNet_Beta(X, y, Adj, lambda0, alpha0,method="APGD",gamma=1000, niter=2000, crit_beta=1e-4, crit_obj=1e-8)
plot(beta_hat_APGD)
```
<img src="man/figures/Beta_hat_APGD.jpeg" width="100%" />

``` r
library("CVXR")
beta_hat_CVX <- HuberNet_Beta(X, y, Adj, lambda0, alpha0,method="CVX")
plot(beta_hat_CVX)
```
<img src="man/figures/Beta_hat_CVX.jpeg" width="100%" />

- **HuberENET**: uber loss function along with Elastic Net penalty function.

``` r
lambda0 = 200
alpha0 = 0.5
beta_hat_APGD <- HuberENET_Beta(X, y, lambda0, alpha0, method="APGD")
library("CVXR")
beta_hat_CVX <- HuberENET_Beta(X, y, lambda0, alpha0,method="CVX")
```


- **HuberLasso**: Huber loss function along with Lasso penalty function.

``` r
lambda0 = 200
beta_hat_APGD <- HuberLasso_Beta(X, y, lambda0, method="APGD")
library("CVXR")
beta_hat_CVX <- HuberLasso_Beta(X, y, lambda0, method="CVX")
```

- **ENET**: Mean square error loss function along with Elastic Net penalty function.

``` r
lambda0 = 200
alpha0 = 0.5
beta_hat_APGD <- ENET_Beta(X, y, lambda0, alpha0, method="APGD")
library("CVXR")
beta_hat_CVX <- ENET_Beta(X, y, lambda0, alpha0,method="CVX")
```

- **Lasso**: Mean square error loss function along with Lasso penalty function.

``` r
lambda0 = 200
beta_hat_APGD <- Lasso_Beta(X, y, lambda0, method="APGD")
library("CVXR")
beta_hat_CVX <- Lasso_Beta(X, y, lambda0, method="CVX")
```

- **Net**: Mean square error loss function along with Network-based penalty function.

``` r
lambda0 = 200
alpha0 = 0.5
beta_hat_APGD <- Net_Beta(X, y, Adj, lambda0, alpha0,method="APGD")
library("CVXR")
beta_hat_CVX <- Net_Beta(X, y, Adj,lambda0, alpha0,method="CVX")
```


### 3. Calculate Selection Probabilities by APGD





<!-- badges: end -->









