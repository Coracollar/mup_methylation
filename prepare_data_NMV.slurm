#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=5
#SBATCH --time=0-30:00:00
#SBATCH --mem=800G

source ~/anaconda3/etc/profile.d/conda.sh

conda activate r_envDSSnew

Rscript prepare_data_NMV.R
