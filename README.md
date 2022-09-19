<!-- README.md is generated from README.Rmd. Please edit that file -->

# TGPred v.0.1.0 (Python Version)

<!-- badges: start -->

Python version of TGPred contains six efficient methods for predicting target genes of a transcription factor by integrating statistics, machine learning, and optimization:

- **HuberNet**: Huber loss function along with Network-based penalty function;
- **HuberLasso**: Huber loss function along with Lasso penalty function;
- **HuberENET**: Huber loss function along with Elastic Net penalty function;
- **MSENet**: Mean square error loss function along with Network-based penalty function;
- **MSELasso**: Mean square error loss function along with Lasso penalty function;
- **MSEENET**: Mean square error loss function along with Elastic Net penalty function;
- **APGD**: The Accelerated Proximal Gradient Descent (APGD) algorithm to solve the above six penalized regression models.

We also have [**R version**](https://github.com/xueweic/TGPred), please see the following link for the guideline of R version https://github.com/xueweic/TGPred.

&emsp;&emsp;


## Installation

***Please use Python Version 3***

**Step 1.** Download 'requrirements.txt' file for installing the requirements packages in your terminal:

``` r
pip install -r requirements.txt
```
&emsp;

**Step 2.** Simple install TGPred package by runing command:

``` r
pip install TGPred
```
&emsp;

**Step 3.** Test your TGPred in python:

``` r
import TGPred
```

If there is no Error, you have installed TGPred package successfully!

&emsp; &emsp;

**Step 4.** Play with function in TGPred referred by [TGPred_Guide](https://github.com/tobefuture/TGPred/blob/main/TGPred_Guide.pdf).

&emsp; &emsp;

## Examples

For quick testing, you can download this files and play with TGPred_example.py.

&emsp; &emsp;

## Reference
Xuewei Cao<sup>+</sup>, Ling Zhang<sup>+</sup>, Mingxia Zhao, Cheng He, Kui Zhang, Sanzhen Liu, Qiuying Sha*, Hairong Wei*. TGPred：Efficient methods for predicting target genes of a transcription factor by integrating statistics, machine learning, and optimization.

<sub> <sup>+</sup> These authors have contributed equally to this work </sub>

**Any questions?** xueweic_AT_mtu_DOT_edu, lingzhan_AT_mtu_DOT_edu

&emsp;&emsp;

## Functions

Please refer [TGPred_Guide](https://github.com/tobefuture/TGPred/blob/main/TGPred_Guide.pdf).
