# Double Robust Variance Estimation

### Bonnie E. Shook-Sa, Paul N. Zivich, Chanhwa Lee, Keyi Xue, Rachael K. Ross, Jessie K. Edwards, Jeffrey S. A. Stringer, and Stephen R. Cole

**Citation**: Shook-Sa BE, Zivich PN, Lee C, Ross RK, Edwards JK, Stringer JSA, Cole SR. "Double Robust Variance Estimation" *arXiv* arXiv:2404.16166.
[![arXiv](https://img.shields.io/badge/arXiv-2404.16166-b31b1b.svg)](https://arxiv.org/abs/2404.16166)
--------------------------------

## Abstract

Doubly robust estimators have gained popularity in the field of causal inference due to their ability to provide consistent point estimates when either an outcome or exposure model is correctly specified. However, the influence function based variance estimator frequently used with doubly robust estimators is only consistent when both the outcome and exposure models are correctly specified. Here, use of M-estimation and the empirical sandwich variance estimator for doubly robust point and variance estimation is demonstrated. Simulation studies illustrate the properties of the influence function based and empirical sandwich variance estimators. Estimators are applied to data from the Improving Pregnancy Outcomes with Progesterone (IPOP) trial to estimate the effect of maternal anemia on birth weight among women with HIV. In the example, birth weights if all women had anemia were estimated to be lower than birth weights if no women had anemia, though estimates were imprecise. Variance estimates were more stable under varying model specifications for the empirical sandwich variance estimator than the influence function based variance estimator. 
--------------------------------

## File Manifesto

### example
The `example/` path contains an example dataset and R and Python code for analyzing the data using the methods presented in the paper.

### sims
The `sims/` path contains R programs to replicate the simulation study presented in the paper.
