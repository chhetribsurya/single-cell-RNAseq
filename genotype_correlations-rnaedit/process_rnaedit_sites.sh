#!/bin/bash

set -x  # Debugging

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
    local outdir="/data/abattle4/surya/datasets/for_lakshmi/results-new/eipl_1_all_auto"
    mkdir -p "$outdir"

    # Extract filename prefix from the VCF file
    local vcf_basename=$(basename "$vcf_file" .vcf.gz)

    echo "Starting process for VCF file: $vcf_file"

    # Function to convert RNA editing sites to BED format
    convert_to_bed() {
        local input_file=$1
        local output_file=$2
        echo "Converting $input_file to BED format..."
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
        echo "Created BED file: $output_file"
    }

    # Function to merge and process BED files
    merge_and_process_bed() {
        local file1=$1
        local file2=$2
        local output_file=$3
        echo "Merging and processing BED files..."
        tail -n+2 "$file1"| cut -f1-3 > "${file1}.filt"
        tail -n+2 "$file2" | cut -f1-3  > "${file2}.filt"
        cat "${file1}.filt" "${file2}.filt" > "${outdir}/${vcf_basename}_merged_rna_edit_sites.filt.bed"
        awk 'BEGIN {OFS="\t"} {if (NR != 1) $1=gensub("chr", "", "g", $1); print}' "${outdir}/${vcf_basename}_merged_rna_edit_sites.filt.bed" > "$output_file"
        echo "Processed merged BED file: $output_file"
    }

    # Function to intersect BED file with VCF
    findoverlaps() {
        local bed_file=$1
        local vcf_file=$2
        local output_file=$3
        echo "Intersecting RNAedit file..."
        bedtools intersect -a "$vcf_file" -b "$bed_file" > "$output_file"
        echo "Created intersected file: $output_file"
    }

    # Convert RNA editing sites to BED format
    convert_to_bed "$darned_rnaedit_sites" "${outdir}/${vcf_basename}_darned.bed"
    convert_to_bed "$redi_rnaedit_sites" "${outdir}/${vcf_basename}_redi.bed"

    # Merge and process BED files
    merge_and_process_bed "${outdir}/${vcf_basename}_darned.bed" "${outdir}/${vcf_basename}_redi.bed" "${outdir}/${vcf_basename}_processed.bed"

    # Remove duplicates
    echo "Removing duplicates from processed BED file..."
    sort -u "${outdir}/${vcf_basename}_processed.bed" > "${outdir}/${vcf_basename}_unique.bed"
    echo "Created unique BED file: ${outdir}/${vcf_basename}_unique.bed"

    # Intersect with VCF file
    findoverlaps "${outdir}/${vcf_basename}_unique.bed" "$vcf_file" "${outdir}/${vcf_basename}_filtered_snps.out.vcf"

    # Count RNA edit sites
    local redit_count=$(wc -l < "${outdir}/${vcf_basename}_filtered_snps.out.vcf")
    echo "Total RNA edit sites present in dataset: $redit_count"

    # Exclude intersections and keep non-overlapping regions
    echo "Creating final filtered VCF file..."
    bedtools intersect -v -a "$vcf_file" -b "${outdir}/${vcf_basename}_unique.bed" -header > "${outdir}/${vcf_basename}_filtered_snps.vcf"
    echo "Created non-overlapping (RNAedit filtered) VCF file: ${outdir}/${vcf_basename}_filtered_snps.vcf"

    # Compress and index the final VCF file
    echo "Compressing and indexing the final VCF file..."
    bcftools view "${outdir}/${vcf_basename}_filtered_snps.vcf" -O z -o "${outdir}/${vcf_basename}_filtered_snps.vcf.gz"
    bcftools index "${outdir}/${vcf_basename}_filtered_snps.vcf.gz"
    echo "Final output file ready: ${outdir}/${vcf_basename}_filtered_snps.vcf.gz"

    # Return the path to the final output file
    echo "${outdir}/${vcf_basename}_filtered_snps.vcf.gz"
}

# Example usage of the function
#vcf_file="/path/to/your_vcf_file.vcf.gz"

gencove_vcf="/data/abattle4/surya/datasets/for_lakshmi/results-new/eipl_1_all_auto/eipl_1_all_modified.vcf.gz"
reference_vcf="/data/abattle4/surya/datasets/for_lakshmi/results-new/eipl_1_all_auto/HPSI0114i-eipl_1.wec.gtarray.HumanCoreExome-12_v1_0.20141111.genotypes_modified.vcf.gz"

processed_reference_vcf=$(process_rna_edit_sites "$reference_vcf")
processed_gencove_vcf=$(process_rna_edit_sites "$gencove_vcf")

echo "Final gencove output vcf file: $processed_gencove_vcf"
echo "Final reference output vcf file: $processed_reference_vcf"
