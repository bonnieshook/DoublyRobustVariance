#!/bin/bash

#SBATCH -p general
#SBATCH -N 1
#SBATCH --mem=1g
#SBATCH -n 1
#SBATCH -t 12:00:00
#SBATCH --mail-type=END,FAIL,REQUEUE,STAGE_OUT
#SBATCH --mail-user=<email address>
module load r/4.1.3

#set number of sims 
N_SIMS=5000

SUBDIR_NAME=N800_SD400
sbatch --output=/dev/null --error=/dev/null --time=2:00:00 --job-name=02_Compile_N800 --wait R CMD BATCH --no-save --no-restore "--args $N_SIMS $SUBDIR_NAME" 02_compile_sims_v2.r 02_Compile_N800.Rout





