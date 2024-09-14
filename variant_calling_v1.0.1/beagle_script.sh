#!/usr/bin/bash

# BEAGLE Function Description:
#   This function automates the process of running the BEAGLE tool for genomic data
#   imputation and phasing across all human chromosomes (1 to 22 and X). It provides
#   flexibility in operation by allowing the user to choose between phasing only,
#   imputation only, or both. Outputs are organized into separate directories based on
#   the operation performed, and each output file is specifically named to reflect
#   the chromosome and the type of operation.

# Parameters:
#   vcf_file: Path to the VCF file containing genomic data. This should be the full
#             path to the VCF file that you want to process.
#   ref_panel: Path to the reference panel VCF file. This should be a VCF file that
#              acts as a reference panel for the imputation process.
#   genetic_map_dir: Directory containing chromosome-specific genetic map files.
#                    This directory should include genetic map files named in the
#                    format 'plink.{chrom}.GRCh37.map', where {chrom} can be any
#                    human chromosome number or 'X'.
#   output_dir: Base directory for storing output files. This directory will contain
#               subdirectories named 'impute_results' or 'phased_results' depending on
#               the operation performed.
#   chr_prefix: Boolean flag (true or false) indicating whether chromosome names in the
#               genetic map files start with 'chr'. Set this to true if the map files
#               use names like 'chr1', 'chr2', etc.
#   impute_option: A string that specifies the operation mode. It can be 'phase' for
#                  phasing only, 'impute' for imputation only, or 'both' for performing
#                  both operations.

# Usage Example:
#   To run the function with imputation set to both phasing and imputation, for a VCF file
#   without chromosome prefixes in the genetic map file names, use the following command:

#run_beagle_for_all_chromosomes "path/to/Zihe_final_output.vcf.gz" \
#                               "path/to/reference_panel.vcf.gz" \
#                               "path/to/genetic_map_dir" \
#                               "path/to/output_directory" \
#                               false \
#                               "both"

# This command processes the specified VCF file for all chromosomes, using the given
# reference panel and genetic maps, performing both phasing and imputation, and stores
# the results in appropriately named directories and files within the specified output directory.

# Function to run BEAGLE for all chromosomes with expanded imputation control
run_beagle_for_all_chromosomes() {
    local vcf_file="$1"
    local ref_panel="$2"
    local genetic_map_dir="$3"
    local output_dir="$4"
    local chr_prefix="$5" # Boolean to include 'chr' prefix in chromosome names
    local impute_option="$6" # Can be phase, impute, both (equals false, true and false/true) 
    echo ${output_dir}

    local impute_actions=()
    case "$impute_option" in
        "phase") impute_actions=("false") ;;
        "impute") impute_actions=("true") ;;
        "both") impute_actions=("true" "false") ;;
    esac

    for impute in "${impute_actions[@]}"; do
        local result_folder="impute_results"
        if [ "$impute" = "false" ]; then
            result_folder="phased_results"
        fi
        local outdir="${output_dir}/${result_folder}"
        mkdir -p "${outdir}"

        # Process all chromosomes (1-22 and X)
        for chromosome in {1..22} X; do
            local chr_name=$chromosome
            [ "$chr_prefix" = true ] && chr_name="chr${chromosome}"

            local map_file="${genetic_map_dir}/plink.${chr_name}.GRCh37.map"
            if [ ! -f "$map_file" ]; then
                echo "Genetic map file not found: $map_file"
                continue  # Skip missing map files
            fi

            local file_suffix="_chr${chr_name}"
            if [ "$impute" = "false" ]; then
                file_suffix="_phased_chr${chr_name}"
            elif [ "$impute" = "true" ]; then
                file_suffix="_imputed_chr${chr_name}"
            fi

            #local file_suffix="_chr${chr_name}${impute:0:1}"
            local vcf_outbase="${outdir}/${vcf_file##*/}${file_suffix}"
            local vcf_outbase="${outdir}/$(basename ${vcf_file} .vcf)${file_suffix}"
            local output_vcf=${vcf_outbase}.vcf.gz
            if [[ ! -s "$output_vcf" ]]; then
                local cmd="beagle \
                    chrom=$chr_name \
                    gt=$vcf_file \
                    ref=$ref_panel \
                    map=$map_file \
                    out=$vcf_outbase \
                    impute=$impute \
                    nthreads=10 \
                    gp=true"

                local action=$( [ "$impute" = "true" ] && echo "Imputation" || echo "Phasing" )
                execute_and_time "BEAGLE $action Chromosome $chr_name" "$action genetic data for chromosome $chr_name using BEAGLE." "$cmd"
                # Index the VCF file with Tabix
                echo "Indexing VCF file with Tabix: $output_vcf"
                tabix_cmd="tabix -p vcf $output_vcf"
                eval "$tabix_cmd"
            fi
        done
    done
}


ref_panel="$REFERENCE_PANEL"
genetic_map_dir="$GENETIC_MAP"
beagle_outdir=$(run_beagle_for_all_chromosomes "$passed_variant_vcf" "$ref_panel" "$genetic_map_dir" "$output_dir" false "both")
echo -e "\n\nBeagle ouput vcf files in dir: $beagle_outdir"




