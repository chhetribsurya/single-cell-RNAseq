#!/bin/bash

# Usage information
usage() {
    echo "Usage: $0 -o <output_directory> -s <'sample1 sample2 ...'>"
    echo "  -o: Output directory for processed files"
    echo "  -s: Space-separated list of sample names (without extension)"
    echo "Fastq files should be named as <sample_name>_1.fastq.gz and <sample_name>_2.fastq.gz"
    exit 1
}


# Parse command line options
while getopts ":o:s:" opt; do
    case ${opt} in
        o)
            output_dir=$OPTARG
            ;;
        s)
            samples=($OPTARG)
            ;;
        \?)
            echo "Invalid Option: -$OPTARG" 1>&2
            usage
            ;;
        :)
            echo "Option -$OPTARG requires an argument." 1>&2
            usage
            ;;
    esac
done


# Check if output directory is set
if [ -z "${output_dir}" ]; then
    echo "Error: Output directory not specified."
    usage
fi

# Check if samples are provided
if [ ${#samples[@]} -eq 0 ]; then
    echo "Error: No samples provided."
    usage
fi


# Create necessary directories
create_dir_if_not_exists() {
    local dir=$1
    if [ ! -d "$dir" ]; then
        echo "Creating directory: $dir"
        mkdir -p "$dir"
    else
        echo "Directory already exists: $dir"
    fi
}

create_dir_if_not_exists "${output_dir}/data"
create_dir_if_not_exists "${output_dir}/output"


# Environment setup
# Define the input directory
FASTQ_DIR="/path_to_fastq_files"

# Define the paths for tools (if direct input)
BWA="/path_to_bwa"
SAMTOOLS="/path_to_samtools"
PICARD="/path_to_picard"
GATK="/path_to_gatk"

# If conda environment is activated. Or given the
# module load (ml) of bwa,samtools, picard, gatk
BWA="bwa"
SAMTOOLS="samtools"
PICARD="picard"
GATK="gatk"

# Define reference files
REFERENCE_GENOME="/path_to_homo_sapiens_assembly38.fasta"
DBSNP="/path_to_dbsnp_146.hg38.vcf.gz"


# Time and execute a command
execute_and_time() {
    local step=$1
    local description=$2
    local command=$3
    local start_time=$(date +%s)

    echo "Starting $step: $description"
    eval $command
    local status=$?

    local end_time=$(date +%s)
    if [ $status -ne 0 ]; then
        echo "Error in $step. Exiting."
        exit $status
    else
        echo "$step completed successfully in $((end_time - start_time)) seconds."
    fi
}


# Function to index reference genome
index_reference_genome() {
    local ref_genome=$1
    execute_and_time "Index Reference Genome" "Indexes the reference genome using BWA." "
        $BWA index ${ref_genome}"
}


# Function for aligning genome with BWA
align_genome() {
    local sample=$1
    local fastq1=$2
    local fastq2=$3
    local bam_output="${output_dir}/output/${sample}.bam"
    execute_and_time "Align Genome" "Aligns raw sequencing data to the human genome using BWA." "
        $BWA mem -M -t 2 \
        ${REFERENCE_GENOME} \
        ${fastq1} \
        ${fastq2} | \
        $SAMTOOLS view -b -h -o ${bam_output} -"
    echo "${bam_output}"
}


# Function for sorting BAM file with Picard
sort_bam() {
    local sample=$1
    local bam_input=$2
    local sorted_bam="${output_dir}/output/${sample}.sort.bam"
    execute_and_time "Sort BAM" "Sorts the BAM file using Picard." "
        $PICARD -Xmx7g SortSam \
        I=${bam_input} \
        O=${sorted_bam} \
        SORT_ORDER=coordinate \
        CREATE_INDEX=True"
    echo "${sorted_bam}"
}


# Function for marking duplicate reads with Picard
mark_duplicates() {
    local sample=$1
    local sorted_bam=$2
    local dedup_bam="${output_dir}/output/${sample}.sort.dup.bam"
    execute_and_time "Mark Duplicates" "Marks duplicate reads in the BAM file using Picard." "
        $PICARD -Xmx7g MarkDuplicates \
        I=${sorted_bam} \
        O=${dedup_bam} \
        M=${output_dir}/output/${sample}_dup_metrics.txt"
    echo "${dedup_bam}"
}


# Function for Base Quality Recalibration with GATK
base_quality_recalibration() {
    local sample=$1
    local dedup_bam=$2
    local recal_bam="${output_dir}/output/${sample}.sort.dup.bqsr.bam"
    execute_and_time "Base Quality Recalibration" "Performs base quality score recalibration with GATK." "
        ${GATK} --java-options \"-Xmx7g\" BaseRecalibrator \
        -I ${dedup_bam} \
        -R ${REFERENCE_GENOME} \
        --known-sites ${DBSNP} \
        -O ${output_dir}/output/${sample}_recal_data.table && \
        ${GATK} --java-options \"-Xmx7g\" ApplyBQSR \
        -I ${dedup_bam} \
        -R ${REFERENCE_GENOME} \
        --bqsr-recal-file ${output_dir}/output/${sample}_recal_data.table \
        -O ${recal_bam}"
    echo "${recal_bam}"
}


# Function for Variant Calling with GATK HaplotypeCaller
variant_calling() {
    local sample=$1
    local recal_bam=$2
    local gvcf="${output_dir}/output/${sample}.g.vcf.gz"
    execute_and_time "Variant Calling" "Calls variants using GATK HaplotypeCaller." "
        ${GATK} --java-options \"-Xmx7g\" HaplotypeCaller \
        -I ${recal_bam} \
        -R ${REFERENCE_GENOME} \
        -ERC GVCF \
        -O ${gvcf}"
    echo "${gvcf}"
}


# Function to combine GVCFs
combine_gvcfs() {
    local gvcf_list=("$@") # Array of GVCF paths
    local combined_gvcf="${output_dir}/output/cohort.g.vcf.gz"
    local gvcf_args=""
    for gvcf in "${gvcf_list[@]}"; do
        gvcf_args+=" -V $gvcf"
    done

    execute_and_time "Combine GVCFs" "Combines multiple GVCF files into a single GVCF with GATK CombineGVCFs." "
        ${GATK} --java-options \"-Xmx7g\" CombineGVCFs \
        -R ${REFERENCE_GENOME} \
        ${gvcf_args} \
        -O ${combined_gvcf}"
    echo "${combined_gvcf}"
}


# Function to perform genotyping
genotyping() {
    local combined_gvcf=$1
    local vcf_output="${output_dir}/output/cohort.vcf.gz"

    execute_and_time "Genotyping" "Performs genotyping on the combined GVCF file with GATK GenotypeGVCFs." "
        ${GATK} --java-options \"-Xmx7g\" GenotypeGVCFs \
        -R ${REFERENCE_GENOME} \
        -V ${combined_gvcf} \
        -O ${vcf_output}"
    echo "${vcf_output}"
}


# Function for Variant Quality Score Recalibration (VQSR) for SNPs
variant_quality_score_recalibration() {
    local vcf_input=$1
    local recal_vcf="${output_dir}/output/cohort_filtered.vcf.gz"

    execute_and_time "Variant Quality Score Recalibration" "Applies VQSR for SNP filtering using GATK." "
        ${GATK} --java-options \"-Xmx7g\" VariantRecalibrator \
        -V ${vcf_input} \
        --resource:hapmap,known=false,training=true,truth=true,prior=15.0 ${REFERENCE_GENOME}/hapmap_3.3.hg38.vcf.gz \
        --resource:omni,known=false,training=true,truth=true,prior=12.0 ${REFERENCE_GENOME}/1000G_omni2.5.hg38.vcf.gz \
        --resource:1000G,known=false,training=true,truth=false,prior=10.0 ${REFERENCE_GENOME}/1000G_phase1.snps.high_confidence.hg38.vcf.gz \
        --resource:dbsnp,known=true,training=false,truth=false,prior=7.0 ${DBSNP} \
        -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
        -mode SNP \
        -O ${output_dir}/output/cohort_snps.recal \
        --tranches-file ${output_dir}/output/cohort_snps.tranches && \
        ${GATK} --java-options \"-Xmx7g\" ApplyVQSR \
        -R ${REFERENCE_GENOME} \
        -V ${vcf_input} \
        --recal-file ${output_dir}/output/cohort_snps.recal \
        --tranches-file ${output_dir}/output/cohort_snps.tranches \
        -mode SNP \
        -O ${recal_vcf}"
    echo "${recal_vcf}"
}


# Function to count variants
count_variants() {
    local vcf_input=$1

    execute_and_time "Count Variants" "Counts the number of variants using GATK." "
        ${GATK} CountVariants \
        -V ${vcf_input}"
}


# Function for additional variant filtering
variant_filtering() {
    local vcf_input=$1
    local filtered_vcf="${output_dir}/output/cohort_filtered_final.vcf.gz"

    execute_and_time "Additional Variant Filtering" "Applies additional variant filters to the VCF file." "
        ${GATK} --java-options \"-Xmx7g\" VariantFiltration \
        -R ${REFERENCE_GENOME} \
        -V ${vcf_input} \
        -O ${filtered_vcf} \
        --filter-name \"DPFilter\" --filter-expression \"DP < 10\""
    echo "${filtered_vcf}"
}


# Function to extract PASS variants
extract_pass_variants() {
    local vcf_input=$1
    local pass_vcf="${output_dir}/output/cohort_final_pass.vcf.gz"

    execute_and_time "Extract PASS Variants" "Extracts variants that passed all filters to create the final VCF file." "
        bcftools view \
        -f PASS \
        ${vcf_input} \
        -Oz -o ${pass_vcf} && \
        tabix -p vcf ${pass_vcf}"
    echo "${pass_vcf}"
}


# Function to convert VCF to TSV for easy viewing
convert_vcf_to_tsv() {
    local vcf_input=$1
    local tsv_output="${output_dir}/output/$(basename "$vcf_input" .vcf.gz)_table.tsv"

    execute_and_time "Convert VCF to TSV" "Converts VCF to TSV format using GATK VariantsToTable for easier data viewing." "
        ${GATK} VariantsToTable \
        -V ${vcf_input} \
        -F CHROM -F POS -F REF -F ALT -F QUAL -GF GT \
        -O ${tsv_output}"
    echo "${tsv_output}"
}


# Function to generate an interactive HTML report
generate_html_report() {
    local sample=$1
    local vcf_input=$2
    local bam_input=$3 # Assuming BAM input is required for the report
    local report_output="${output_dir}/output/${sample}_report.html"

    execute_and_time "Generate HTML Report" "Generates an interactive HTML report for variant visualization." "
        jigv --sample ${sample} \
        --sites ${vcf_input} \
        --fasta ${REFERENCE_GENOME} \
        ${bam_input} > ${report_output}"
    echo "${report_output}"
}

# Function to generate QC reports
generate_qc_reports() {
    local sample=$1
    local fastq1=$2
    local fastq2=$3
    local qc_output_dir="${output_dir}/output/QC_reports/${sample}"

    create_dir_if_not_exists "$qc_output_dir"
    execute_and_time "Generate QC Reports" "Generates QC reports for raw fastq data using FastQC and aggregates them using MultiQC." "
        fastqc -o ${qc_output_dir} ${fastq1} ${fastq2} &&
        multiqc -o ${qc_output_dir} ${qc_output_dir}"
}

# Function to archive and compress results
archive_results() {
    local sample=$1
    local archive_file="${output_dir}/output/${sample}_pipeline_results.tar.gz"

    execute_and_time "Archive and Compress Results" "Archives and compresses the output files for easy storage and transfer." "
        tar -czvf ${archive_file} -C ${output_dir}/output ."
    echo "${archive_file}"
}


# Index the reference genome
index_reference_genome "${REFERENCE_GENOME}"


# Pipeline Execution
for sample_name in "${samples[@]}"; do

    echo "Processing sample: $sample_name"

    # Define paths to input FastQ files
    fastq_file_1="${FASTQ_DIR}/${sample_name}_1.fastq.gz"
    fastq_file_2="${FASTQ_DIR}/${sample_name}_2.fastq.gz"

    # Check if FastQ files exist
    if [ ! -f "$fastq_file_1" ] || [ ! -f "$fastq_file_2" ]; then
        echo "Fastq files for sample $sample_name not found. Skipping."
        continue
    fi

    # Align the genome with fastqs
    bam_output=$(align_genome "$sample_name" "$fastq_file_1" "$fastq_file_2")

    # Sort BAM file
    sorted_bam=$(sort_bam "$sample_name" "$bam_output")

    # Mark duplicates
    dedup_bam=$(mark_duplicates "$sample_name" "$sorted_bam")

    # Base quality recalibration
    recal_bam=$(base_quality_recalibration "$sample_name" "$dedup_bam")

    # Variant calling
    gvcf=$(variant_calling "$sample_name" "$recal_bam")

    # Variant counts before additional filters
    count_variants "${gvcf}"

done


# Cohort-wide processing
# Find all GVCF files in the output directory
cohort_gvcfs=($(find "${output_dir}/output" -name "*.g.vcf.gz"))

# Check if there are multiple GVCF files for cohort-wide analysis
if [ ${#cohort_gvcfs[@]} -gt 1 ]; then

    # Combine GVCFs
    combined_gvcf=$(combine_gvcfs "${cohort_gvcfs[@]}")

    # Genotyping
    vcf_output=$(genotyping "$combined_gvcf")

    # Variant Quality Score Recalibration (VQSR)
    filtered_vcf=$(variant_quality_score_recalibration "$vcf_output")

    # Variant counts
    count_variants "${filtered_vcf}"

    # Additional Variant filtering
    final_vcf=$(variant_filtering "$filtered_vcf")

    # Extract PASS Variants
    pass_vcf=$(extract_pass_variants "$final_vcf")

    # Convert VCF to TSV for easier viewing
    tsv_output=$(convert_vcf_to_tsv "$pass_vcf")

    # Generate HTML report
    report_output=$(generate_html_report "$sample_name" "$pass_vcf" "$recal_bam")

elif [ ${#cohort_gvcfs[@]} -eq 1 ]; then

    echo "Single GVCF file found. Proceeding with single-sample analysis."
    # Continue with single-sample analysis steps ...

else

    echo "No GVCF files found in the output directory. Skipping cohort-wide steps."

fi


# Generate QC reports and archive results for each sample
for sample_name in "${samples[@]}"; do

    # Define paths to input FastQ files
    fastq_file_1="${FASTQ_DIR}/${sample_name}_1.fastq.gz"
    fastq_file_2="${FASTQ_DIR}/${sample_name}_2.fastq.gz"

    # Generate QC reports
    generate_qc_reports "$sample_name" "$fastq_file_1" "$fastq_file_2"
    archive_file=$(archive_results "$sample_name")

done

echo "Pipeline completed for all samples."

