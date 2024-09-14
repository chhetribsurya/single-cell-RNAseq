#!/bin/bash

#SBATCH --job-name=genotype_correlation
#SBATCH --output=%LOG_DIR%/%B_PARAM%_%A_%a.log
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=10
#SBATCH --mem=35G
#SBATCH --partition=shared

# Module load from server
module load anaconda

# Call the main script with dynamically generated args
#./genotype_correlation_main.sh \
./genotype_correlation_mainplus.sh \
    -e gatk_env \
    -o %OUTPUT_DIR% \
    -b %B_PARAM% \
    -g %G_PARAM% \
    -r %R_PARAM%
