#!/bin/bash

# SLURM Script for SeqTK Automation
# This script submits a job to the SLURM queue to automate seqtk operations on specified input files.

#SBATCH --job-name=seqtk_automate          # Job name
#SBATCH --output=seqtk_automate_%j.log     # Standard output and error log (with job ID)
#SBATCH --time=14:00:00                    # Time limit hrs:min:sec
#SBATCH --cpus-per-task=16                 # Number of CPU cores per task
#SBATCH --mem=62G                          # Total memory limit
#SBATCH --partition=defq                   # Partition to submit to, replace with your partition

# Execute the seqtk automation script
./seqtk_automate.sh
