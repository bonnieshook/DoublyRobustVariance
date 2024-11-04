# Double Robust Variance Estimation with Parametric Working Models

### Bonnie E. Shook-Sa, Paul N. Zivich, Chanhwa Lee, Keyi Xue, Rachael K. Ross, Jessie K. Edwards, Jeffrey S. A. Stringer, and Stephen R. Cole

**Citation**: Shook-Sa BE, Zivich PN, Lee C, Ross RK, Edwards JK, Stringer JSA, Cole SR. "Double Robust Variance Estimation with Parametric Working Models" *arXiv* arXiv:2404.16166.
[![arXiv](https://img.shields.io/badge/arXiv-2404.16166-b31b1b.svg)](https://arxiv.org/abs/2404.16166)
--------------------------------

## Abstract

Doubly robust estimators have gained popularity in the field of causal inference due to their ability to provide consistent point estimates when either an outcome or exposure model is correctly specified. However, for nonrandomized exposures the influence function based variance estimator frequently used with doubly robust estimators of the average causal effect is only consistent when both working models (i.e., outcome and exposure models) are correctly specified. Here, the empirical sandwich variance estimator and the nonparametric bootstrap are demonstrated to be doubly robust variance estimators. That is, they are expected to provide valid estimates of the variance leading to nominal confidence interval coverage when only one working model is correctly specified. Simulation studies illustrate the properties of the influence function based, empirical sandwich, and nonparametric bootstrap variance estimators in the setting where parametric working models are assumed. Estimators are applied to data from the Improving Pregnancy Outcomes with Progesterone (IPOP) study to estimate the effect of maternal anemia on birth weight among women with HIV.
--------------------------------

## File Manifesto

### example
The `example/` path contains an example dataset and R and Python code for analyzing the data using the methods presented in the paper.

### sims
The `sims/` path contains R programs to replicate the simulation study presented in the paper.
