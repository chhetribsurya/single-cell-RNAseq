#!/bin/bash

# Check if the correct number of arguments was provided
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <input_vcf.gz> <output_vcf.gz> <chromosomes>"
    exit 1
fi

# Assign command line arguments to variables
INPUT_VCF=$1
OUTPUT_VCF=$2
CHROMOSOMES=$3  # Comma-separated list of chromosomes to include

# Step 1: Filter out specified chromosomes
echo "Filtering chromosomes: $CHROMOSOMES..."
bcftools view -r $CHROMOSOMES -Oz -o intermediate_output.vcf.gz $INPUT_VCF

# Step 2: Extract the header from the filtered file
echo "Extracting header..."
bcftools view -h intermediate_output.vcf.gz > original_header.txt

# Step 3: Modify the header to remove specific contig lines for Y and MT
echo "Modifying header..."
grep -v -E '##contig=<ID=(Y|MT|NC_007605|hs37d5|GL[0-9]+\.?[0-9]*),' original_header.txt > modified_header.txt

#grep -v -E '##contig=<ID=(Y|MT|NC_007605|hs37d5|GL[0-9]+),' original_header.txt > modified_header.txt

#grep -v -E '##contig=<ID=(Y|MT|GL)' original_header.txt > modified_header.txt
#grep -v -E '##contig=<ID=(Y|MT|GL|NC_007605|hs37d5),' original_header.txt > modified_header.txt

#grep -v -E '##contig=<ID=(Y|MT),' original_header.txt > modified_header.txt


# Step 4: Apply the new header to the filtered VCF
echo "Applying new header..."
bcftools reheader -h modified_header.txt -o temp_output.vcf.gz intermediate_output.vcf.gz

# Step 5: Compress and index the final VCF file using tabix
echo "Compressing and indexing final VCF file..."
bcftools view -Oz temp_output.vcf.gz > $OUTPUT_VCF
tabix -p vcf $OUTPUT_VCF

# Optional: Clean up intermediate files
echo "Cleaning up..."
rm original_header.txt modified_header.txt temp_output.vcf.gz intermediate_output.vcf.gz

echo "Process completed successfully."

