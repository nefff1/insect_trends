#!/bin/bash

#SBATCH --job-name=insect_trends

#SBATCH --time=255:00:00
#SBATCH -A node
#SBATCH -p node
#SBATCH --qos=normal
#SBATCH --output=/home/neff/out_bin/bin_insect_trends_%A_%a.txt
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=2000
#SBATCH --array=1-72 # example for dragonflies

SCRIPT=/home/neff/R_Occupancy_detection_models.R

srun  singularity /home/neff/ubuntu_container Rscript --vanilla $SCRIPT $SLURM_ARRAY_TASK_ID 2000




