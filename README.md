
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
Xuewei Cao<sup>+</sup>, Ling Zhang<sup>+</sup>, Kui Zhang, Sanzhen Liu, Qiuying Sha<sup>*</sup>, Hairong Wei<sup>*</sup>. HuberNet function for interfering target genes of regulatory genes using high-throughput gene expression data.

<sub><sup> <sup>+</sup> These authors have contributed equally to this work </sup></sub>

**Any questions**? xueweic_AT_mtu_DOT_edu, lingzhan_AT_mtu_DOT_edu

## Example

### Simulated data

``` r
library(APGD)
```


<img src="man/figures/README-pressure-1.png" width="100%" />


