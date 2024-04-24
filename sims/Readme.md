Double Robust Variance Estimation

R code for simulation study in Section 4 of the manuscript.
These programs implement the simulation study described in Shook-Sa et al, "Double Robust Variance Estimation" in R

The simulations are implemented in the following order:

1. Calculate the true ACE empirically using 00_DRsim_Truth_SD400.R
1. Generate and analyze data for the simulation study using the 01_DRsim_one_iter_SD400.R program, called with the 01_Run_DR_Sims_N800_SD400.sh shell script. This process was repeated for each of the other sample size and sigma scenarios presented in the Supporting Information.
   Note that this program calls the 00_Functions.R program, which contains functions for computing each estimator of interest along with estimated standard errors.
2. Combine the results from all iterations of the simulation using the 02_compile_sims_v2.r program, called with the shell script 02_Compile_Sims_all.sh. This program creates a final dataset with one row for each iteration of the simulations.
3. Calculate summary measures using 03_combine_sim_tables800_SD400.R and generate the variance ratio figure using 04_var_figure_v2_800_SD400.R.

Developed by: Bonnie Shook-Sa
