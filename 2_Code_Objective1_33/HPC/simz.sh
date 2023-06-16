#!/bin/bash
#SBATCH --export=ALL
#SBATCH --ntasks=1                               # Number of cores
#SBATCH --nodes=1                                # Number of nodes for the cores
#SBATCH --job-name="simz.sh"                     # Name of the script
#SBATCH --time=0-08:00                           # Runtime in D-HH:MM format i
#SBATCH --partition=cpu                          # Partition to submit to (cpu or gpu)
#SBATCH --mem=2                                  # Memory pool for all CPUs i
#SBATCH --mail-type=ALL                          # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=j.munozavila@umcutrecht.nl   # Email to which notifications will be sent
#SBATCH --array=1-10                            # Array job=start-end, or =x,y,z, where the job number is assigned to $SLURM_ARRAY_TASK_ID

# File to which standard out will be written
#SBATCH --output=/home/julius_te/jmunoz/Results_z1/%j.out

# File to which standard err will be written
#SBATCH --error=/home/julius_te/jmunoz/Results_z1/%j.err

# Path to R and the local package library
PATH=/hpc/local/CentOS7/julius_te/R-4.0.4/bin/:$PATH
R_libs=/hpc/local/CentOS7/julius_te/R-4.0.4/lib64/R/library


# check which R version we use
which R

# Rscript makes R run a given script, with the following arguments:
# --vanilla specifies the default options for Rscript
# path to user made R script
# parameters passed to the user made R script.

Rscript --vanilla /home/julius_te/jmunoz/Run_zika/run_zika1.R $SLURM_ARRAY_TASK_ID

# end script


