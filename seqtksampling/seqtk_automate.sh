#!/bin/bash

# Seqtk Automation Script
# This script performs seqtk sampling on specified donor fastq files and compresses the output.
# Usage: ./seqtk_automate.sh

# Ensure seqtk is installed and available in PATH.
eval "$(conda shell.bash hook)"
conda activate seqtk

# Directories for input and output
input="/data/abattle4/lakshmi/cuomo_2020/pipeline_files/concat_files"
output="/data/abattle4/surya/datasets/for_lakshmi/seqtk_outputs"
#output="/data/abattle4/lakshmi/cuomo_2020/pipeline_files/seqtk_outputs"
#mkdir -p $output


# Array of donors
donors=("xojn_3" "letw_1" "sojd_3" "poih_4" "joxm_1")

# Loop through each donor to generate seqtk outputs
for donor in "${donors[@]}"; do
    echo "Processing donor: ${donor}"

    # Define input files for R1 and R2
    input_r1="${input}/${donor}_R1.fastq.gz"
    input_r2="${input}/${donor}_R2.fastq.gz"

    # Check if input files exist
    if [[ ! -f "$input_r1" ]] || [[ ! -f "$input_r2" ]]; then
        echo "Input files for ${donor} not found. Skipping."
        continue
    fi

    # Seqtk sub-sample command
    echo "Running seqtk on ${donor}"
    seqtk sample -s 42 "$input_r1" 80000000 > "${output}/${donor}_R1_finalc.fastq"
    seqtk sample -s 42 "$input_r2" 80000000 > "${output}/${donor}_R2_finalc.fastq"
done

# Gzip all created files
echo "Compressing output files..."
for fastq in ${output}/*.fastq; do
    echo "Processing fastq: $fastq"
    gzip -f "$fastq"
done

echo "Seqtk processing completed."
