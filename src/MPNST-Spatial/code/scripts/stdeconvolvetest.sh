#!/usr/bin/env bash
#SBATCH -p cpu
#SBATCH -n 30
#SBATCH --mem=700G
#SBATCH --gres=gpu:0
#SBATCH -t 72:0:00


source /etc/profile.d/modules.sh

module load R

Rscript /home/mli110/MPNST-Spatial/scripts/scttestscript.R











