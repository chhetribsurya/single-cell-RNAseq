#!/bin/bash

#SBATCH --job-name=genotype_correlation
#SBATCH --output=%LOG_DIR%/%B_PARAM%_%A_%a.log
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=6
#SBATCH --mem=20G
#SBATCH --partition=shared

# Module load from server
# module load anaconda

# Call the main script with dynamically generated args
#./genotype_correlation_main.sh \
echo -e "\nCommand: ./genotype_correlation_mainplus.sh \
    -e /home/schhetr1/anaconda3/envs/bcftools \
    -o %OUTPUT_DIR% \
    -b %B_PARAM% \
    -g %G_PARAM% \
    -r %R_PARAM% \n"

./genotype_correlation_mainplus.sh \
    -e /home/schhetr1/anaconda3/envs/bcftools \
    -o %OUTPUT_DIR% \
    -b %B_PARAM% \
    -g %G_PARAM% \
    -r %R_PARAM%
