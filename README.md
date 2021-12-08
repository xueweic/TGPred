
<!-- README.md is generated from README.Rmd. Please edit that file -->

# APGD v.0.1.0

<!-- badges: start -->
<!-- badges: end -->

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




<img src="man/figures/README-pressure-1.png" width="100%" />


