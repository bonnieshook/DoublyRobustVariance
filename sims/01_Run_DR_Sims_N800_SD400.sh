#!/bin/bash

#SBATCH -p general
#SBATCH -N 1
#SBATCH --mem=2g
#SBATCH -n 1
#SBATCH -t 12:00:00
#SBATCH --mail-type=END,FAIL,REQUEUE,STAGE_OUT
#SBATCH --mail-user=<email address>
module load r/4.1.3

#define variables
N_SIMS=5000 
SUBDIR_NAME=N800_SD400
N=800
BOOTS=1000
SIGMA=400

sbatch --output=/dev/null --error=/dev/null --time=2:00:00 --array=1-$N_SIMS --job-name=01_Run_DR_Sims_N800 --wait R CMD BATCH --no-save --no-restore "--args $SUBDIR_NAME $N $BOOTS $SIGMA" 01_DRsim_one_iter_SD400.R 01_DRsim_one_iter_SD400.Rout






