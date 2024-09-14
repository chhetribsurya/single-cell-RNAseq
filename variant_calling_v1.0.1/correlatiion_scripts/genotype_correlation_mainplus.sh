#!/usr/bin/env bash

# --------------------------------------------------------------------------------
# Script Name: genotype_correlation_main.sh
# Pipeline Contact: Surya Chhetri/Ashton Omdahl/Lakshmi Chanemougam
# Description:
#   This script processes VCF genomic data files to compute genotype correlations.
#   It utilizes tools such as bcftools and PLINK to filter and analyze VCF files.
#
# Usage:
#   ./genotype_correlation_main.sh [options]
#
# Notes:
#   - Ensure all dependencies are installed and accessible in your environment.
#   - This script is part of the single cell eQTL project.
#
# --------------------------------------------------------------------------------

# Detailed Usage:
# ./script_name.sh [-e ENV_NAME] [-o OUTPUT_DIR] [-b BASE_OUTPUT_FILE] [-g GENCOVE_FILE] [-r REFERENCE_FILE]
# Options:
# -e, --env            Set the conda environment name (default: sam_bcfenv)
# -o, --output-dir     Set the base directory for output files
# -b, --base-output-file Set the base name for output files (default: output_file)
# -g, --gencove-file   Set the Gencove output VCF filename
# -r, --reference-file Set the Reference genotype VCF filename


# Function to display help message
function show_help() {
    echo "Usage: $0 [-e ENV_NAME] [-o OUTPUT_DIR] [-b BASE_OUTPUT_FILE] [-g GENCOVE_FILE] [-r REFERENCE_FILE]"
    echo "Options:"
    echo "  -e, --env              Set the conda environment name (default: sam_bcfenv)"
    echo "  -o, --output-dir       Set the base directory for output files"
    echo "  -b, --base-output-file Set the base name for output files (default: output_file)"
    echo "  -g, --gencove-file     Set the Gencove output VCF filename"
    echo "  -r, --reference-file   Set the Reference genotype VCF filename"
    echo "  -h, --help             Display this help and exit"
}

# # Default parameters
env_name="sam_bcfenv"
env_name="bcftools"

base_output_file="outputfile"
# base_output_file="zihe_1"

# Function to parse command line arguments
function parse_arguments() {
    while [ "$#" -gt 0 ]; do
        case $1 in
            -e|--env)
                env_name="$2"; shift 2;;
            -o|--output-dir)
                output_dir="$2"; shift 2;;
            -b|--base-output-file)
                base_output_file="$2"; shift 2;;
            -g|--gencove-file)
                gencove_file="$2"; shift 2;;
            -r|--reference-file)
                reference_file="$2"; shift 2;;
            -h|--help)
                show_help
                exit 0;;
            *)
                echo "Unknown option: $1" >&2
                show_help
                exit 1;;
        esac
    done
}

# Parse the command-line arguments
parse_arguments "$@"

# Initialize Conda
eval "$(conda shell.bash hook)"
conda activate "$env_name"

# Ensure output directory exists
output_dir="$output_dir/$base_output_file"
mkdir -p "$output_dir"

# Record time log function
function record_log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
}

# Function to modify the ID column of VCF files
function set_vcf_id() {
    local vcf_file="$1"
    local basename_vcf="$(basename ${vcf_file%.vcf.gz}_modified.vcf.gz)"
    local modified_vcf="${output_dir}/${basename_vcf}"
    bcftools annotate --set-id '%CHROM:%POS' "$vcf_file" -O z -o "$modified_vcf"
    echo "$modified_vcf"
}

# Define full paths for files
gencove_vcf="$gencove_file"
reference_vcf="$reference_file"

# Define the function for VCF processing
process_vcf() {
    local output_dir="$1"
    local base_output_file="$2"
    local gencove_vcf="$3"
    local reference_vcf="$4"

    # Processing steps with debug logs
    record_log "Starting VCF processing."
    
    # Modify VCF ID for Gencove and Reference VCF files
    modified_gencove_vcf=$(set_vcf_id "$gencove_vcf")
    modified_reference_vcf=$(set_vcf_id "$reference_vcf")

    # Define full paths for modified files
    gencove_vcf="$modified_gencove_vcf"
    reference_vcf="$modified_reference_vcf"

    # Step 1: Filter VCF of Gencove output by 'PASS'
    filtered_gencove_vcf="${output_dir}/${base_output_file}_gencove_filtered.vcf.gz"
    record_log "Filtering Gencove VCF by 'PASS'."
    bcftools view -f "PASS" -m2 -M2 --types snps --output-type z \
    --output "$filtered_gencove_vcf" "$gencove_vcf" && record_log "Gencove VCF filtered successfully at $(date '+%Y-%m-%d %H:%M:%S')."

    # Step 2: Extract bi-allelic SNPs with score cutoff of 0.15
    filtered_reference_vcf="${output_dir}/${base_output_file}_ref_filtered.vcf.gz"
    record_log "Extracting bi-allelic SNPs with score cutoff of 0.15."
    bcftools view -m2 -M2 --types snps -i 'GC>0.15' --output-type z \
    --output "$filtered_reference_vcf" "$reference_vcf" && record_log "Bi-allelic SNPs extracted successfully at $(date '+%Y-%m-%d %H:%M:%S')."

    # Step 3: Extract SNP ids and reference allele from filtered VCF
    ref_alleles_tsv="${output_dir}/${base_output_file}_ref_filtered_alleles.tsv"
    record_log "Extracting SNP ids and reference allele from filtered VCF."
    bcftools query -f '%ID\t%REF\n' -i 'BAF>=0.01' \
    "$filtered_reference_vcf" > "$ref_alleles_tsv" && record_log "SNP ids and reference allele extracted successfully at $(date '+%Y-%m-%d %H:%M:%S')."

    # Step 4: Load appropriate version of plink and make bed files for filtered Gencove SNPs
    ml plink/2.00a2.3
    # ml plink/1.90b6.4

    duplicate_ids="${output_dir}/duplicateIDs.txt"
    bcftools query -f '%ID\n' "$filtered_gencove_vcf" | sort | uniq -d > "$duplicate_ids"

    harmonized_bed="${output_dir}/${base_output_file}_gencove_HARMONIZED"
    record_log "Making bed files for filtered Gencove SNPs."

    plink2 --silent --vcf "$filtered_gencove_vcf" \
    --exclude "$duplicate_ids" \
    --ref-allele "$ref_alleles_tsv" \
    --make-bed --out "$harmonized_bed" && record_log "Bed files for filtered Gencove SNPs created successfully at $(date '+%Y-%m-%d %H:%M:%S')."

    # # If goal is to generate .pgen, .psam, and .pvar files instead, 
    # # you would use the --pfile flag rather than --make-bed

    # plink2 --vcf "$filtered_gencove_vcf" \
    #        --exclude "$duplicate_ids" \
    #        --ref-allele "$ref_alleles_tsv" \
    #        --pfile "$harmonized_bed" \
    #        --output-chr MT

    # Step 5: Extract just the snps that match based on SNP IDs
    ref_rsid_tsv="${output_dir}/${base_output_file}_hipsci_reference_RSIDS.tsv"
    record_log "Extracting matching SNP IDs."
    cut -f 1 "$ref_alleles_tsv" > "$ref_rsid_tsv" && record_log "Matching SNP IDs extracted successfully at $(date '+%Y-%m-%d %H:%M:%S')."

    # Step 6: Extract the matching SNP IDs with Plink
    harmonized_subset_vcf="${output_dir}/${base_output_file}_gencove_harmonized_subsetted"
    record_log "Extracting the matching SNP IDs with Plink."

    plink2 --silent --bfile "$harmonized_bed" \
           --extract "$ref_rsid_tsv" \
           --recode vcf \
           --out "$harmonized_subset_vcf" && record_log "Matching SNP IDs extracted successfully with Plink at $(date '+%Y-%m-%d %H:%M:%S')."

    # Check if plink2 command was successful
    harmonized_subset_vcf="${output_dir}/${base_output_file}_gencove_harmonized_subsetted.vcf"

    if [ $? -eq 0 ]; then
        echo "VCF file created successfully."
        # Gzip the VCF file
        gzip "$harmonized_subset_vcf"
        echo "Gzipped VCF file created at: ${harmonized_subset_vcf}.gz"
    else
        echo "Error occurred in plink2 command."
        exit 1
    fi

    # Step 7: Extract the positions of SNPs from reference
    filtered_positions_tsv="${output_dir}/${base_output_file}_ref_filtered_positions.tsv"

    bcftools query -f '%CHROM\t%POS\n' \
    --output "$filtered_positions_tsv" \
    "$filtered_reference_vcf"

    if [ $? -eq 0 ]; then
        # echo "Positions of SNPs extracted successfully at $(date '+%Y-%m-%d %H:%M:%S')."
        # Optional: Call to record_log function
        record_log "Positions of SNPs extracted successfully at $(date '+%Y-%m-%d %H:%M:%S')."
    else
        echo "Error occurred during bcftools query."
        exit 1
    fi

    # Step 8: Quality control - Indexing the Gencove SNPs PASS filter output
    record_log "Indexing the Gencove SNPs PASS filter output."
    bcftools index "$filtered_gencove_vcf" && record_log "Indexing completed at $(date '+%Y-%m-%d %H:%M:%S')."

    # Step 9: Quality control - Extract and check data
    sanity_check_file="${output_dir}/${base_output_file}_sanity_check"
    record_log "Performing sanity check."
    bcftools query --regions-file "$filtered_positions_tsv" \
    -i 'AF>=0.01' -f '%CHROM\tPOS\t%INFO\t%FILTER[\t%DS]\n' \
    "$filtered_gencove_vcf" > "$sanity_check_file" && record_log "Sanity check completed at $(date '+%Y-%m-%d %H:%M:%S')."

    # Check for "LOWCONF" entries in the file
    if grep -q "LOWCONF" "$sanity_check_file"; then
        # If "LOWCONF" is found, log the message
        record_log "Low confidence entries found in sanity check at $(date '+%Y-%m-%d %H:%M:%S')."
    else
        # If "LOWCONF" is not found, output the opposite message
        echo "No low confidence entries found."
    fi

    # Step 10: Extract dosages from final VCF files
    record_log "Extracting dosages from final VCF files."
    bcftools query -f '%CHROM\t%POS[\t%GT]\n' "$filtered_reference_vcf" > "${output_dir}/${base_output_file}_ref_filtered.csv"
    bcftools query -f '%CHROM\t%POS[\t%GT]\n' "${harmonized_subset_vcf}.gz" > "${output_dir}/${base_output_file}_gencove_harmonized_subsetted.csv" && record_log "Dosages extracted successfully at $(date '+%Y-%m-%d %H:%M:%S')."

    record_log "VCF processing completed."

}

# Call the VCF processing function with arguments
process_vcf "$output_dir" "$base_output_file" "$gencove_vcf" "$reference_vcf"
echo -e "\nVCF function call completed ...\n"

# Initialize python conda env or module load python
env="py3_newenv"
eval "$(conda shell.bash hook)"
conda activate $env

# Define paths of the output files from the bash script
ref_output_file="${output_dir}/${base_output_file}_ref_filtered.csv"
imputed_output_file="${output_dir}/${base_output_file}_gencove_harmonized_subsetted.csv"

# Define Python and script path
echo -e "\nPython function call for correlation compute ...\n"
python_script_path="./genotype_correlation_final.py"  # Replace with the actual path to the Python script

# Run the Python script with output files as arguments
python "$python_script_path" "$ref_output_file" "$imputed_output_file" "${base_output_file}"
echo -e "\nCorrelation compute done ...\n"
