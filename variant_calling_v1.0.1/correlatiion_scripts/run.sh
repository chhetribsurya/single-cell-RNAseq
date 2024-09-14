#!/bin/bash

# Read prefixes from the file
mapfile -t prefixes < file_prefixes.txt

# Suffix to be appended to each prefix
suffix="_gatk_pass"

# Output directory
output_dir="/data/abattle4/surya/datasets/for_lakshmi/variant_calling/correlation_genotype_vcfs/results-gatk"

# Base directories for g_param and r_param
#g_base_dir="/data/abattle4/surya/datasets/for_lakshmi/variant_calling/correlation_genotype_vcfs"
#g_base_dir="/data/abattle4/surya/datasets/for_lakshmi/variant_calling/correlation_genotype_vcfs/PASS_VCFs"
g_base_dir="/data/abattle4/surya/datasets/for_lakshmi/variant_calling/correlation_genotype_vcfs/STAR_VERSION/PASS_VCFs"
r_base_dir="/data/abattle4/lakshmi/cuomo_2020/genotypes"

# Path to the SLURM script template
slurm_template="./genotype_sbatch_submit_job.sh"

# Check if g_base_dir exists, if not, exit with error message
if [ ! -d "$g_base_dir" ]; then
    echo "Error: Directory $g_base_dir does not exist. Please check the directory path."
    exit 1
fi

# Check if r_base_dir exists, if not, exit with error message
if [ ! -d "$r_base_dir" ]; then
    echo "Error: Directory $r_base_dir does not exist. Please check the directory path."
    exit 1
fi

# Check if output directory exists
if [ ! -d "$output_dir" ]; then
    mkdir -p "$output_dir"
fi

# Log files directory
log_dir="./log_files"

# Create log files directory if it does not exist
if [ ! -d "$log_dir" ]; then
    mkdir -p "$log_dir"
fi

# Get the length of the prefixes array
length=${#prefixes[@]}

# Submit jobs
for ((i=0; i<length; i++)); do
    
    # Find the g_param and r_param in the directories
    b_prefix=${prefixes[i]}
    b_param=${b_prefix}${suffix}
    g_param=$(find $g_base_dir -type f -name "*${b_prefix}*.vcf.gz" -print -quit)
    r_param=$(find $r_base_dir -type f -name "*${b_prefix}*.genotypes.vcf.gz" -print -quit)

    # Create a temp SLURM script from the template
    temp_slurm_script="slurm_$b_param.sh"

    cat $slurm_template | \
        sed "s|%B_PARAM%|$b_param|g" | \
        sed "s|%G_PARAM%|$g_param|g" | \
        sed "s|%R_PARAM%|$r_param|g" | \
        sed "s|%OUTPUT_DIR%|$output_dir|g" | \
        sed "s|%LOG_DIR%|$log_dir|g" > $temp_slurm_script

    # Submit the temporary SLURM script
    sbatch $temp_slurm_script

    # Optionally remove the temporary script after submission
    rm $temp_slurm_script

    #sbatch --export=g_param="$g_param",r_param="$r_param",output_dir="$output_dir",log_dir="$log_dir" genotype_sbatch_submit_jobs.sh
 
done
