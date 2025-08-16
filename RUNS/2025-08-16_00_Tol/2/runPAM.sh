#!/bin/bash

#SBATCH -p general
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --cpus-per-task=1
#SBATCH --mem=64g
#SBATCH -t 07-00:00:00

module load matlab/2024b
matlab -nodesktop -nosplash -singleCompThread -r runner_postprocessing -logfile PAM.out


