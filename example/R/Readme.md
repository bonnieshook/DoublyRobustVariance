Double Robust Variance Estimation with Parametric Working Models

R code for applying each of the three doubly robust estimators with three variance estimation approaches to the data in exampledata.csv

The dataset 'exampledata.csv' contains the following variables:

part.id = participant id (1-n)
Z1 = continuous covariate
Z2 = binary covariate
Z3 = binary covariate
X = binary exposure
Y = continuous outcome 

The program 01_Analysis.R demonstrates how to analyze the example dataset using the three doubly robust estimators demonstrated in the paper. 
Correct model specification is based on the data generating mechanism described in Section 4 of the manuscript.
This program calls the 00_Functions.R program, which contains the functions for each estimator. 

Developed by: Bonnie Shook-Sa (Thanks to Bradley Saul for streamlining the weighted regression AIPW functions)