#!/usr/bin/env bash
#SBATCH -p cpu
#SBATCH -n 30
#SBATCH --mem=500G
#SBATCH --gres=gpu:0
#SBATCH -t 100:0:00
#SBATCH --mail-user=mli110@jh.edu
#SBATCH --mail-type=ALL


source /etc/profile.d/modules.sh

module load R

Rscript /home/mli110/MPNST-Spatial/scripts/nmf_all_samples_k6to20.R