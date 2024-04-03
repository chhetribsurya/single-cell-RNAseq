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

process_rna_edit_sites() {

    # Module loads
    ml bedtools/2.30.0
    ml bcftools/1.15.1

    # Input VCF file
    local vcf_file=$1

    # Paths to RNA editing sites files
    local darned_rnaedit_sites="/data/abattle4/surya/datasets/for_lakshmi/data/rnaedit_sites/darned_hg19.txt"
    local redi_rnaedit_sites="/data/abattle4/surya/datasets/for_lakshmi/data/rnaedit_sites/rediPortal_TABLE1_hg19.txt"

    # Output directory
    #local outdir="/data/abattle4/surya/datasets/for_lakshmi/results-new/eipl_1_all_auto"
    local outdir="$output_dir"
    mkdir -p "$outdir"

    # Extract filename prefix from the VCF file
    local vcf_basename=$(basename "$vcf_file" .vcf.gz)

    #echo "Starting process for VCF file: $vcf_file"

    # Function to convert RNA editing sites to BED format
    convert_to_bed() {
        local input_file=$1
        local output_file=$2
        #echo "Converting $input_file to BED format..."
        awk 'BEGIN {OFS="\t"} {
            if (NR==1) {
                $2="start";
                print;
            } else {
                printf "%s\t%s\t%s", $1, $2-1, $2;
                for (i=4; i<=NF; i++) {
                    printf "\t%s", $i;
                }
                printf "\n";
            }
        }' "$input_file" > "$output_file"
        #echo "Created BED file: $output_file"
    }

    # Function to merge and process BED files
    merge_and_process_bed() {
        local file1=$1
        local file2=$2
        local output_file=$3
        #echo "Merging and processing BED files..."
        tail -n+2 "$file1"| cut -f1-3 > "${file1}.filt"
        tail -n+2 "$file2" | cut -f1-3  > "${file2}.filt"
        cat "${file1}.filt" "${file2}.filt" > "${outdir}/${vcf_basename}_merged_rna_edit_sites.filt.bed"
        awk 'BEGIN {OFS="\t"} {if (NR != 1) $1=gensub("chr", "", "g", $1); print}' "${outdir}/${vcf_basename}_merged_rna_edit_sites.filt.bed" > "$output_file"
        #echo "Processed merged BED file: $output_file"
    }

    # Find overlaps
    findoverlaps() {
        local bed_file=$1
        local vcf_file=$2
        local output_file=$3
        #echo "Intersecting RNAedit file..."
        bedtools intersect -a "$vcf_file" -b "$bed_file" > "$output_file"
        #echo "Created intersected file: $output_file"
    }

    # Convert RNA editing sites to BED format
    convert_to_bed "$darned_rnaedit_sites" "${outdir}/${vcf_basename}_darned.bed"
    convert_to_bed "$redi_rnaedit_sites" "${outdir}/${vcf_basename}_redi.bed"

    # Merge and process BED files
    merge_and_process_bed "${outdir}/${vcf_basename}_darned.bed" "${outdir}/${vcf_basename}_redi.bed" "${outdir}/${vcf_basename}_processed.bed"

    # Remove duplicates
    #echo "Removing duplicates from processed BED file..."
    sort -u "${outdir}/${vcf_basename}_processed.bed" > "${outdir}/${vcf_basename}_unique.bed"
    #echo "Merged unique RNA edit file: ${outdir}/${vcf_basename}_unique.bed" >&2

    # Remove CN.* to align bed format and intersect
    cnvfilt_vcf="${outdir}/${vcf_basename}_cnvfiltered.forbedint.vcf.gz"
    bcftools view "$vcf_file" -e 'ALT ~ "<CN.*>"' -Oz -o "$cnvfilt_vcf" 
    bcftools index "$cnvfilt_vcf"

    # Intersect with VCF file
    vcf_file=$cnvfilt_vcf
    findoverlaps "${outdir}/${vcf_basename}_unique.bed" "$vcf_file" "${outdir}/${vcf_basename}_filtered_snps.out.vcf"

    # Count RNA edit sites
    #local redit_count=$(wc -l < "${outdir}/${vcf_basename}_filtered_snps.out.vcf")
    #echo -e "\nRNA editing sites in ${vcf_file}: $redit_count \n" >&2

    # Exclude intersections and keep non-overlapping regions
    #echo "Creating final filtered VCF file..."
    #bedtools intersect -v -a "$vcf_file" -b "${outdir}/${vcf_basename}_unique.bed" -header > "${outdir}/${vcf_basename}_filtered_snps.vcf"
    bedtools intersect -a "$vcf_file" -b "${outdir}/${vcf_basename}_unique.bed" -header > "${outdir}/${vcf_basename}_filtered_snps.vcf"
    echo -e "\nRNA editing sites filtered out VCF file: ${outdir}/${vcf_basename}_filtered_snps.vcf" >&2

    # Compress and index the final VCF file
    #echo "Compressing and indexing the final VCF file..."
    bcftools view "${outdir}/${vcf_basename}_filtered_snps.vcf" -O z -o "${outdir}/${vcf_basename}_filtered_snps.vcf.gz"
    bcftools index "${outdir}/${vcf_basename}_filtered_snps.vcf.gz"
    #echo "Final output file ready: ${outdir}/${vcf_basename}_filtered_snps.vcf.gz"

    # Return the path to the final output file
    echo "${outdir}/${vcf_basename}_filtered_snps.vcf.gz"
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
    
    # Processing and filtering out RNA editing sites
    record_log "Filtering out RNA editing sites  and processing VCFs"
    rnaedit_processed_gencove_vcf=$(process_rna_edit_sites "$gencove_vcf")
    rnaedit_processed_reference_vcf=$(process_rna_edit_sites "$reference_vcf")
    record_log "RNA editing sites filtered out successfully at $(date '+%Y-%m-%d %H:%M:%S')."
    record_log "GENCOVE VCF RNAedit sites removed: $rnaedit_processed_gencove_vcf"
    record_log "REFERENCE VCF RNAedit sites removed: $rnaedit_processed_reference_vcf"

    # Modify VCF ID for Gencove and Reference VCF files
    record_log "Setting VCF IDs"
    modified_gencove_vcf=$(set_vcf_id "$rnaedit_processed_gencove_vcf")
    modified_reference_vcf=$(set_vcf_id "$rnaedit_processed_reference_vcf")
    record_log "VCF IDs set  successfully at $(date '+%Y-%m-%d %H:%M:%S')."
    # modified_gencove_vcf=$(set_vcf_id "$gencove_vcf")
    # modified_reference_vcf=$(set_vcf_id "$reference_vcf")

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
eval "$(conda shell.bash hook)"
#conda activate my-env
conda activate /home/schhetr1/anaconda3/envs/py3-env

# Define paths of the output files from the bash script
ref_output_file="${output_dir}/${base_output_file}_ref_filtered.csv"
imputed_output_file="${output_dir}/${base_output_file}_gencove_harmonized_subsetted.csv"

# Define Python and script path
echo -e "\nPython function call for correlation compute ...\n"
python_script_path="./genotype_correlation_final.py"  # Replace with the actual path to the Python script

# Run the Python script with output files as arguments
python "$python_script_path" "$ref_output_file" "$imputed_output_file" "${base_output_file}"
echo -e "\nCorrelation compute done ...\n"
