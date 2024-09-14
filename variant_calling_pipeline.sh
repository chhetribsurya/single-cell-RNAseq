#!/bin/bash

# Script to extract a specific chromosome from a genome FASTA file,
# build a reference index, and align FASTQ files to this reference, producing a BAM file.

# Requires: samtools, bowtie2, bwa

# Usage:
# chmod +x align_to_chromosome.sh
# ./align_to_chromosome.sh --genome-fasta <genome_fasta> --chromosome <chromosome> --output-dir <output_dir> --fastq1 <fastq1> [--fastq2 <fastq2>] [--aligner <bowtie2|bwa>]

# Display detailed usage instructions
usage() {
    echo "Usage: $0 --genome-fasta <genome_fasta> --chromosome <chromosome> --output-dir <output_dir> --fastq1 <fastq1> [--fastq2 <fastq2>] [--aligner <bowtie2|bwa>]"
    echo ""
    echo "Arguments:"
    echo " --genome-fasta <genome_fasta> Path to the genome FASTA file. This file should be in FASTA format."
    echo " --chromosome <chromosome> Specific chromosome or region to be extracted (e.g., 'chr20', 'chrX:1000-5000', "chr1, chr2, chr3" or 'all')."
    echo " --output-dir <output_dir> Directory where output files will be saved. Defaults to the current directory."
    echo " --fastq1 <fastq1> Path to the first FASTQ file. Required for both single-end and paired-end reads."
    echo " --fastq2 <fastq2> Path to the second FASTQ file. Optional, only for paired-end reads."
    echo " --aligner <bowtie2|bwa> Aligner to use. Options: 'bowtie2' or 'bwa'. Defaults to 'bowtie2'."
    echo " --sample <sample> Name of the sample being processed. Defaults to 'outputfile'."
    echo " --genome-index Enable genome indexing mode only and exit immediately after completion. No value needed."
    echo ""
    echo ""
    echo " --jumpToBam <true/false> If specified, the script will process bam files directly,"
    echo "                          skipping the alignment from FastQ files. Use this flag if"
    echo "                          you already have bamfiles as input instead of FastQ/CRAM files."
    echo " --bam_file <bam file> Name of the sample being processed. Defaults to 'outputfile'."
    echo ""
    echo "--jumpToCram <true/false> If specified, the script will process CRAM files directly,"
    echo "                           skipping the alignment from FastQ files. Use this flag if"
    echo "                          you already have CRAM files as input instead of FastQ/BAM files."
    echo " --cram_file <cram file> Name of the sample being processed. Defaults to 'outputfile'."
    echo ""
    echo "Notes:"
    echo "  - This script assumes that all necessary tools (like BWA, SAMtools, GATK) are"
    echo "    installed and properly configured in the execution environment."
    echo ""
    echo "Example with paired-end fastqs:"
    echo " $0 --genome-fasta path/to/genome.fa --chromosome chr1 --output-dir path/to/output --fastq1 path/to/read1.fastq --fastq2 path/to/read2.fastq --aligner bwa"
    exit 1
}


# Default values for flags
jumpToCram=false
jumpToBam=false

# Default parameters
output_dir="."
aligner="bowtie2"
sample="resultOutfile"


# Parse flag-based command-line arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
        --genome-fasta) genome_fasta="$2"; shift 2 ;;
        --chromosome) chromosome="$2"; shift 2 ;;
        --output-dir) output_dir="$2"; shift 2 ;;
        --fastq1) fastq1="$2"; shift 2 ;;
        --fastq2) fastq2="$2"; shift 2 ;;
        --aligner) aligner="$2"; shift 2 ;;
        --sample) sample="$2"; shift 2 ;;
        --genome-index) genome_index="true"; shift ;;
        -h|--help) usage ;;
        *) echo "Unknown option: $1"; usage ;;
    esac
done


# Create main output dir if doesn't exist
mkdir -p "$output_dir/samples"

# Checks if log file exists to avoid log appends
mkdir -p "$output_dir/pipelogs"
logfile="$output_dir/pipelogs/pipe-${sample}.log"
if [ -f "$logfile" ]; then
    rm "$logfile"
fi

echo -e "\nLOG file for sample $sample at $logfile"

# Logging function
log() {
    echo -e "[$(date +'%Y-%m-%d %H:%M:%S')] $@" >> "$logfile"
}


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

#create_dir_if_not_exists "${output_dir}/output"


# Define reference files
DATA_DIR="/scratch16/abattle4/surya/datasets/for_gatk/resource_bundles/GRCh37"
DBSNP="${DATA_DIR}/dbsnp_138.b37.vcf.gz"
HAPMAP="${DATA_DIR}/hapmap_3.3.b37.vcf.gz"
OMNI="${DATA_DIR}/1000G_omni2.5.b37.vcf.gz"
ONEKGENOME="${DATA_DIR}/1000G_phase1.snps.high_confidence.b37.vcf.gz"
#REFERENCE_GENOME="${DATA_DIR}/hs37d5.fa"
#CONTIGS_CONVERSION_FILE="${DATA_DIR}/contigs_conversion.txt"
GENETIC_MAP="${DATA_DIR}/plink_map"
REFERENCE_PANEL="${DATA_DIR}/ALL.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"

# Command execution
cmd_execution() {
    local step="$1"
    local description="$2"
    local command="$3"
    local start_time=$(date +%s)
    
    log "Executing $step: $description"
    log "Command: $command"
    eval "$command"
    local status=$?
    local end_time=$(date +%s)
    local duration=$((end_time - start_time))
    
    if [ $status -ne 0 ]; then
        log "Error in $step. Status: $status. Exiting."
        exit $status
    else
        log "$step completed successfully in $duration seconds.\n"
    fi
}


# Execute a command, log its activity, and time the execution
execute_and_time() {
    local step="$1"
    local description="$2"
    local command="$3"
    local max_attempts=1  # Maximum attempts including the first try
    local attempt=0
    local success=0

    log "Executing $step: $description"

    while [ $attempt -lt $max_attempts ]; do
        ((attempt++))
        local start_time=$(date +%s)
        log "Attempt $attempt: Command: $command"
        eval "$command"
        local status=$?
        local end_time=$(date +%s)
        local duration=$((end_time - start_time))

        if [ $status -ne 0 ]; then
            log "Error in $step during attempt $attempt. Status: $status. Duration: $duration seconds."
            if [ $attempt -eq $max_attempts ]; then
                log "Failed after $attempt attempts. Exiting."
                exit $status
            fi
            log "Retrying $step..."
        else
            log "$step completed successfully in $duration seconds on attempt $attempt.\n"
            success=1
            break
        fi
    done

    if [ $success -ne 1 ]; then
        log "Unable to complete $step successfully after $max_attempts attempts."
        exit 1
    fi
}


# Function to index the genome with samtools faidx
index_genome() {
    local genome_fasta="$1"
    local output_dir="$2"
    #local outdir="${output_dir}/$sample"
    local output_fasta="${output_dir}/fullchrom.fa"
    echo ${output_fasta}
    if [[ ! -s "${output_fasta}.fai" ]]; then
        local cmd1="cp $genome_fasta $output_fasta"
        local cmd2="samtools faidx $output_fasta"
        execute_and_time "Copy Genome" "Copying genome fasta file." "$cmd1"
        execute_and_time "Index Genome" "Indexing the genome fasta file." "$cmd2"
        #eval "$cmd2" && execute_and_time "Index Genome" "Indexing the genome fasta file." "$cmd2"
    fi
}


# Function to extract chromosome using samtools faidx
extract_chromosome() {
    local genome_fasta="$1"
    local chromosome="$2"
    local output_dir="$3"
    #local outdir="${output_dir}/$sample"
    local strname=$(echo "$chromosome" | tr ' ' '_')
    local output_fasta="${output_dir}/chromset.${strname}.fa"
    echo "$output_fasta"
    if [[ ! -s "$output_fasta" ]]; then
        local cmd1="samtools faidx $genome_fasta $chromosome > $output_fasta"
        local cmd2="samtools faidx $output_fasta"
        execute_and_time "Extract Chromosome" "Extracting $chromosome from the genome." "$cmd1"
        execute_and_time "Index Extracted Chromosome" "Indexing the extracted chromosome fasta file." "$cmd2"
        #execute_and_time "Extract Chromosome" "Extracting $chromosome from the genome." "$cmd1" && eval "$cmd2"
    fi
}


build_reference_index() {
    local ref_genome="$1"
    local aligner="$2"

    if [[ "$aligner" == "bwa" ]]; then
        # Make sure the command is enclosed correctly in quotes
        local index_prefix="${ref_genome}"
        echo "$index_prefix"  # For consistency in output
        if [[ ! -s "${index_prefix}.bwt" ]]; then
            local cmd="bwa index $ref_genome"
            execute_and_time "Build Reference Index" "Building reference index using BWA." "$cmd"
        fi

    elif [[ "$aligner" == "bowtie2" ]]; then
        local index_prefix="${ref_genome%.fa}_bowtie2"
        echo "$index_prefix"  # For consistency in output
        if [[ ! -s "$index_prefix" ]]; then
            local cmd="bowtie2-build $ref_genome $index_prefix"
            execute_and_time "Build Reference Index" "Building reference index using Bowtie2." "$cmd"
        fi
    else
        echo "Error: Unsupported aligner specified."
        exit 1
    fi
}


# Function for processing CRAM files: convert to BAM, sort, and index
process_cram() {
    local sample=$1
    local cram_file="${output_dir}/data/${sample}.cram"
    local sorted_bam="${output_dir}/output/${sample}_sorted.bam"

    execute_and_time "Process CRAM" "Converts CRAM to BAM, sorts, and indexes the BAM file." "
        $SAMTOOLS view -T $REFERENCE_GENOME -b ${cram_file} | \
        $SAMTOOLS sort -o ${sorted_bam} - &&
        $SAMTOOLS index ${sorted_bam}"
    echo "${sorted_bam}"
}


# Function to parse FASTQ to extract RG IDs
parse_fastq_headers() {
    local file="$1"
    local first_line=$(zcat "$file" | head -n 1)

    # Assuming standard Illumina header format: @Instrument:RunID:FlowCellID:Lane:Tile:X:Y ReadNum:FilterControl:Index
    local instrument=$(echo "$first_line" | cut -d ':' -f1 | cut -d '@' -f2| xargs| tr " " "_")
    local run_id=$(echo "$first_line" | cut -d ':' -f2)
    local flow_cell_id=$(echo "$first_line" | cut -d ':' -f3)
    local lane=$(echo "$first_line" | cut -d ':' -f4)
    #local tile=$(echo "$first_line" | cut -d ':' -f5 | sed 's/\/[12]$//')  # Remove the /1 or /2 from the tile information
    #local sample_id=$(basename "$file" | sed 's/_R[12].*//')  # Extract sample ID from filename

    # Construct unique ID
    local unique_id="${instrument}_${run_id}_${flow_cell_id}_${lane}"
    echo "$unique_id"

    # Construct the RG line
    #local rg_line="@RG\tID:${sample_id}_${unique_id}\tSM:$sample_id\tLB:$instrument\tPL:ILLUMINA"
}


# Function to align reads
align_reads() {
    local ref_index="$1"
    local fastq1="$2"
    local fastq2="$3"
    local aligner="$4"
    local output_dir="$5"
    local rgid="$6"
    local sam_output="${output_dir}/$sample/${sample}.sam"
    local bam_output="${output_dir}/$sample/${sample}.bam"
    echo "$bam_output"

    if [[ "$aligner" == "bwa" ]]; then

        #cmd="bwa mem -M -t 10 -R "@RG\tID:ERR2970877.1_HS31_25475_3_1305_9226_81781\tSM:zihe_1\tLB:ERR2970877.1\tPL:ILLUMINA" $ref_index $fastq1 ${fastq2:+$fastq2} | samtools view -bS - > $bam_output"
        cmd="bwa mem -M -t 10 -R \"@RG\\tID:${sample}_${rgid}\\tSM:${sample}\\tLB:${rgid}\\tPL:ILLUMINA\" $ref_index $fastq1 ${fastq2:+$fastq2} | samtools view -bS - > $bam_output"

    elif [[ "$aligner" == "bowtie2" ]]; then

        if [[ -z "$fastq2" ]]; then
            cmd="bowtie2 -x $ref_index -U $fastq1 | samtools view -bS - > $bam_output"
        else
            cmd="bowtie2 -x $ref_index -1 $fastq1 -2 $fastq2 | samtools view -bS - > $bam_output"
        fi

    fi

    if [[ ! -s "$bam_output" ]]; then
        execute_and_time "Align Reads" "Aligning reads to the reference using $aligner." "$cmd"
    fi
}


# Function for sorting BAM file with Picard
sort_bam() {
    local sample="$1"
    local bam_input="$2"
    local output_dir="$3"
    local sorted_bam="${output_dir}/${sample}/${sample}.sort.bam"
    echo "${sorted_bam}"
    if [[ ! -s "${sorted_bam}" ]]; then
        local cmd="picard -Xmx50g SortSam \
            I=${bam_input} \
            MAX_RECORDS_IN_RAM=3000000 \
            O=${sorted_bam} \
            SORT_ORDER=coordinate \
            CREATE_INDEX=True"
        
        execute_and_time "Sort BAM" "Sorts the BAM file using Picard." "$cmd"
    fi
}


# Function for marking duplicate reads with Picard
mark_duplicates() {
    local sample="$1"
    local sorted_bam="$2"
    local output_dir="$3"
    local dedup_bam="${output_dir}/${sample}/${sample}.sort.dup.bam"
    echo "${dedup_bam}"
    if [[ ! -s "$dedup_bam" ]]; then
        local cmd="picard -Xmx50g MarkDuplicates \
            I=${sorted_bam} \
            MAX_RECORDS_IN_RAM=3000000 \
            O=${dedup_bam} \
            M=${output_dir}/${sample}_dup_metrics.txt"
         execute_and_time "Mark Duplicates" "Marks duplicate reads in the BAM file using Picard." "$cmd"
    fi
}


# Function to create a FASTA dictionary file using GATK
create_fasta_dictionary() {
    #local reference_fasta="${REFERENCE_GENOME}"
    local reference_fasta="$1"
    local outdict_file="${reference_fasta%.fa}.dict"
    #if [ -f "${outdict_file}" ]; then
    #    rm ${outdict_file}
    #fi
    echo "$outdict_file"
    if [[ ! -s "$outdict_file" ]]; then
        local cmd="gatk CreateSequenceDictionary -R $reference_fasta -O $outdict_file"
        execute_and_time "Create FASTA Dictionary" "Creating a dictionary for the reference genome $reference_fasta." "$cmd"
    fi
}


# Function to index a VCF file using GATK IndexFeatureFile
index_vcf_file() {
    #local vcf_file="${DBSNP}"
    local vcf_file="$1"
    echo "$vcf_file"
    if [[ ! -s "$vcf_file" ]]; then
        local cmd="gatk IndexFeatureFile -I $vcf_file"
        execute_and_time "Index VCF File" "Indexing VCF file $vcf_file." "$cmd"
    fi
}


# Function to re-compress and index a VCF file
recompress_and_index_vcf() {
    local vcf_path="$1"
    local vcf_file=$(basename "$vcf_path")
    local vcf_dir=$(dirname "$vcf_path")

    if [[ ! -s "${vcf_path}.tbi" ]]; then
        # Decompress if already gzipped but not bgzipped
        if file --mime-type "$vcf_path" | grep -q gzip; then
            gunzip "$vcf_path"
        fi

        ## Recompress with bgzip
        bgzip "${vcf_path%.gz}"

        # Index the VCF file
        local bgzipped_vcf="$vcf_dir/${vcf_file%.gz}.gz"
        echo "Indexing $bgzipped_vcf"
        local cmd="gatk IndexFeatureFile -I $bgzipped_vcf"
        execute_and_time "Index VCF File" "Indexing recompressed bgzipped VCF file $vcf_file." "$cmd"
        echo "Index created for: $bgzipped_vcf"
    fi
}


# Function to check contigs and conditionally reheader a VCF
check_and_conditionally_reheader_vcf() {
    local reference_fasta="$1"
    local input_vcf="$2"
    local conversion_file="$3"
    local output_vcf="${input_vcf%.vcf.gz}_reheadered.vcf.gz"

    # Extract contigs from reference FASTA and VCF
    local ref_contigs=$(samtools faidx "$reference_fasta" && cut -f1 "$reference_fasta.fai")
    local vcf_contigs=$(bcftools view -h "$input_vcf" | grep "^##contig=<ID=" | sed 's/.*ID=\([^,]*\).*/\1/')

    # Compare contigs
    local differences=$(echo -e "$ref_contigs\n$vcf_contigs" | sort | uniq -u)
    if [ -n "$differences" ]; then
        echo "Discrepancies found in VCF-contigs. Proceeding with reheadering..."
        local cmd="bcftools reheader -f $conversion_file $input_vcf -o $output_vcf"
        execute_and_time "Reheadering VCF" "After VCF-contigs incompatibility reheaders VCF ${input_vcf}." "$cmd"

        local cmd="tabix -p vcf $output_vcf"
        execute_and_time "Indexing Reheadered VCF" "After VCF-contigs incompatibility indexes reheadered VCF ${input_vcf}." "$cmd"

        echo "Reheadered VCF file created at: $output_vcf"
        DBSNP="$output_vcf" # Update/reset the global DBSNP var to the new VCF file

    else
        echo "No discrepancies found in VCF-contigs. Using the original VCF file."
        DBSNP="$input_vcf" # Retain the DBSNP global var to the same original VCF file
    fi
}


# Extract VCF contigs, compare them, and conditionally rename chromosomes
check_and_conditionally_rename_vcf() {
    local reference_fasta="$1"
    local input_vcf="$2"
    local conversion_remapfile="$3"
    local output_vcf="${input_vcf%.vcf.gz}_reheadered.vcf.gz"

    # Extract contigs from VCF using bcftools and process them
    local ref_contigs=$(samtools faidx "$reference_fasta" && cut -f1 "$reference_fasta.fai" | sort | uniq| head -n 10)
    local vcf_contigs=$(bcftools view -h "$input_vcf" | grep "^##contig=<ID=" | sed 's/.*ID=\([^,]*\).*/\1/' | sort | uniq| head -n 10)

    # Count the entries with "chr" prefix
    local ref_has_chr=$(echo "$ref_contigs" | grep -c '^chr')
    local vcf_has_chr=$(echo "$vcf_contigs" | grep -c '^chr')

    # Compare the extracted VCF contigs with reference contigs
    #local differences=$(echo -e "$ref_contigs\n$vcf_contigs" | sort | uniq -u)
    #if [ -n "$differences" ]; then

    # Compare the counts of "chr" prefixes
    if [[ $ref_has_chr -ne $vcf_has_chr ]]; then
        echo "Discrepancy detected in VCF contigs based on 'chr' prefix diff. Proceeding with renaming..."
        local cmd="bcftools annotate --rename-chrs $conversion_remapfile -Oz -o $output_vcf $input_vcf"
        execute_and_time "Rename Chromosomes in VCF" "Renaming chromosomes in VCF file to match reference genome contigs." "$cmd"
        
        local index_cmd="tabix -p vcf $output_vcf"
        execute_and_time "Index VCF" "Indexing the renamed VCF file." "$index_cmd"

        echo "Renamed and indexed VCF file created at: $output_vcf"
        DBSNP="$output_vcf"  # Update the DBSNP to the new VCF file
    else
        echo "No discrepancies found. Using the original VCF file."
        #cp "$input_vcf" "$output_vcf"
        #tabix -p vcf "$output_vcf" # Ensure the copied file is indexed
        DBSNP="$input_vcf"  # Keep DBSNP pointing to the original VCF
    fi
}


# Function for Base Quality Recalibration with GATK
base_quality_recalibration() {
    local sample="$1"
    local dedup_bam="$2"
    local output_dir="$3"
    local outdir="${output_dir}/${sample}/results"
    local recal_bam="${outdir}/${sample}.sort.dup.bqsr.bam"
    echo "${recal_bam}"
    mkdir -p $outdir
    if [[ ! -s "$recal_bam" ]]; then
        local cmd="gatk  BaseRecalibrator \
            -I ${dedup_bam} \
            -R ${REFERENCE_GENOME} \
            --known-sites ${DBSNP} \
            -O ${outdir}/${sample}_recal_data.table && \
            gatk ApplyBQSR \
            -I ${dedup_bam} \
            -R ${REFERENCE_GENOME} \
            --bqsr-recal-file ${outdir}/${sample}_recal_data.table \
            -O ${recal_bam}"
        
        execute_and_time "Base Quality Recalibration" "Performs base quality score recalibration with GATK for sample ${sample}." "$cmd"
    fi
}


post_base_quality_recalibration() {
    local sample="$1"
    local recal_bam="$2"  # This should be the recalibrated BAM file from ApplyBQSR
    local output_dir="$3"
    local outdir="${output_dir}/${sample}/results"
    local recal_table="${outdir}/${sample}_recal_data.table"  # First recalibration table
    local after_recal_table="${outdir}/${sample}_after_recal_data.table"  # New recalibration table to be generated
    local plots="${outdir}/${sample}_recal_plots.pdf"
    mkdir -p "${outdir}"
    if [[ ! -s "${after_recal_table}" ]]; then
        local cmd="gatk BaseRecalibrator \
            -R ${REFERENCE_GENOME} \
            -I ${recal_bam} \
            --known-sites ${DBSNP} \
            -O ${after_recal_table} && \
            gatk AnalyzeCovariates \
            -before ${recal_table} \
            -after ${after_recal_table} \
            -plots ${plots}"

        execute_and_time "Base Quality Recalibration Post-Analysis" "Generating recalibration report and analyzing covariates for sample ${sample}." "$cmd"
        #echo "Recalibration analysis and plots available at: ${plots}"
    fi
}


# Function for Base Quality Recalibration Post-Analysis
post_base_quality_recalibration_only() {
    local sample="$1"
    local dedup_bam="$2"
    local output_dir="$3"
    local outdir="${output_dir}/${sample}/results"
    local recal_table="${outdir}/${sample}_recal_data.table"
    local after_recal_table="${outdir}/${sample}_after_recal_data.table"
    #echo "${after_recal_table}"
    mkdir -p "${outdir}"
    local cmd="gatk BaseRecalibrator \
        -R ${REFERENCE_GENOME} \
        -I ${dedup_bam} \
        --known-sites ${DBSNP} \
        -BQSR ${recal_table} \
        -o ${after_recal_table}"
        
        #-knownSites ${KNOWN_SITES_1} \
        #-knownSites ${KNOWN_SITES_2} \
        #-BQSR ${recal_table} \
        #-o ${after_recal_table}"

    execute_and_time "Base Quality Recalibration Post-Analysis" "Calculating the quality scores after recalibration for sample ${sample}." "$cmd"
}


# Function for Analyzing Covariates
analyze_covariates() {
    local sample="$1"
    local output_dir="$2"
    local outdir="${output_dir}/${sample}/results"
    local recal_table="${outdir}/${sample}_recal_data.table"
    local after_recal_table="${outdir}/${sample}_after_recal_data.table"
    local plots="${outdir}/${sample}_recal_plots.pdf"
    mkdir -p "${outdir}"

    local cmd="gatk AnalyzeCovariates \
        -before ${recal_table} \
        -after ${after_recal_table} \
        -plots ${plots}"

     #-R ${REFERENCE_GENOME} \
    execute_and_time "Analyze Covariates" "Generating plots to analyze the effects of recalibration for sample ${sample}." "$cmd"
    echo "${plots}"
}


# Function for Variant Calling with GATK HaplotypeCaller
variant_calling() {
    local sample="$1"
    local recal_bam="$2"
    local output_dir="$3"  # Add the output directory as a parameter
    local outdir="${output_dir}/${sample}/results"
    local gvcf="${outdir}/${sample}.g.vcf.gz"
    mkdir -p "${outdir}"  # Ensure the output directory exists
    echo "${gvcf}"
    if [[ ! -s "${gvcf}" ]]; then
        local cmd="gatk HaplotypeCaller \
            -I ${recal_bam} \
            -R ${REFERENCE_GENOME} \
            -ERC GVCF \
            --standard-min-confidence-threshold-for-calling 30.0 \
            -O ${gvcf}"

        #-ERC BP_RESOLUTION \
        #-ERC GVCF \
        #-stand-call-conf 30.0 \
        execute_and_time "Variant Calling" "Calls variants using GATK HaplotypeCaller for sample ${sample}." "$cmd"
    fi
}


# Function for genotyping using GATK GenotypeGVCFs, handling both single and multi-sample modes
genotype_gvcf_withnolist() {
    local sample="$1"
    local gvcf="$2"  # Single gVCF for single mode; ignored in multi mode
    local output_dir="$3"
    local mode="$4"  # "single" or "multi"
    local outdir="${output_dir}/${sample}/results"
    local final_vcf="${outdir}/${sample}_final_output.SINGLE.vcf.gz"
    echo "${final_vcf}"
    mkdir -p "${outdir}"  # Ensure the output directory exists

    if [[ ! -s "${final_vcf}" ]]; then
        if [ "${mode}" = "single" ]; then
            #echo "${final_vcf}"
            #echo "Running single-sample genotyping..."
            local cmd="gatk GenotypeGVCFs \
                -R ${REFERENCE_GENOME} \
                -V ${gvcf} \
                -O ${final_vcf}"
            execute_and_time "Single-sample Genotype GVCF" "Genotyping single GVCF for sample ${sample} using GATK GenotypeGVCFs." "$cmd"
        else
            #echo "${final_vcf}"
            #echo "Running multi-sample genotyping..."
            # Assume all gVCF files are in ${outdir} for multi-sample mode
            local combined_gvcf="${outdir}/combined.g.vcf.gz"
            local gvcf_args=$(find ${outdir} -name '*.g.vcf.gz' -exec echo --variant {} \; | tr '\n' ' ')
            local cmd_combine="gatk CombineGVCFs \
                -R ${REFERENCE_GENOME} \
                ${gvcf_args} \
                -O ${combined_gvcf}"
            execute_and_time "Combine GVCFs" "Combining GVCFs for multiple samples using GATK CombineGVCFs." "$cmd_combine"

            local cmd_genotype="gatk GenotypeGVCFs \
                -R ${REFERENCE_GENOME} \
                -V ${combined_gvcf} \
                -O ${final_vcf}"
            execute_and_time "Multi-sample Genotype GVCF" "Genotyping combined GVCF for multiple samples using GATK GenotypeGVCFs." "$cmd_genotype"
        fi
    fi

    #echo "Final VCF generated at: ${final_vcf}"
}


# Function for genotyping using GATK GenotypeGVCFs, handling both single and multi-sample modes
genotype_gvcf() {
    local sample="$1"
    local gvcf="$2"  # Single gVCF for single mode; ignored in multi mode
    local output_dir="$3"
    local mode="$4"  # "single" or "multi"
    local outdir="${output_dir}/${sample}/results"
    local final_vcf="${outdir}/${sample}_final_output.vcf.gz"
    echo "${final_vcf}"
    mkdir -p "${outdir}"  # Ensure the output directory exists

    if [[ ! -s "${final_vcf}" ]]; then
        if [ "${mode}" = "single" ]; then
            #echo "Running single-sample genotyping..."
            local cmd="gatk GenotypeGVCFs \
                -R ${REFERENCE_GENOME} \
                -V ${gvcf} \
                -O ${final_vcf}"
            execute_and_time "Single-sample Genotype GVCF" "Genotyping single GVCF for sample ${sample} using GATK GenotypeGVCFs." "$cmd"
        else
            # Assume all gVCF files are in ${outdir} for multi-sample mode
            local result_dir="$output_dir/multisample_results"
            mkdir -p $result_dir
            local combined_gvcf="${result_dir}/combined.g.vcf.gz"
            if [[ ! -s "${combined_gvcf}" ]]; then
                # Generate a file containing the list of variants
                local variant_list_file="${result_dir}/gvcf_list.txt"
                find "${output_dir}" -name '*.g.vcf.gz' > "${variant_list_file}"
                local cmd_combine="gatk CombineGVCFs \
                    -R ${REFERENCE_GENOME} \
                    --variant ${variant_list_file} \
                    -O ${combined_gvcf}"
                execute_and_time "Combine GVCFs" "Combining GVCFs for multiple samples using GATK CombineGVCFs." "$cmd_combine"

                local cmd_genotype="gatk GenotypeGVCFs \
                    -R ${REFERENCE_GENOME} \
                    -V ${combined_gvcf} \
                    -O ${final_vcf}"
                execute_and_time "Multi-sample Genotype GVCF" "Genotyping combined GVCF for multiple samples using GATK GenotypeGVCFs." "$cmd_genotype"
            fi
        fi
    fi
}


# Function for Variant Quality Score Recalibration (VQSR) for SNPs
variant_quality_score_recalibration() {
    local sample="$1"
    local vcf_input="$2"
    local output_dir="$3"
    local filt_thresh="99"
    local outdir="${output_dir}/${sample}/results"
    local recal_vcf="${outdir}/${sample}_VQSR.snp.recal.${filt_thresh}.vcf"
    local recal_file="${outdir}/${sample}_VQSR.snp.recal"
    local tranches_file="${outdir}/${sample}_VQSR.snp.tranches"
    local rscript_file="${outdir}/${sample}_VQSR.snp.plots.R"
    echo "${recal_vcf}"
    mkdir -p "${outdir}"  # Ensure the output directory exists
    
    if [[ ! -s "${recal_vcf}" ]]; then
        #echo "Output VCF will be: ${recal_vcf}"
        #echo "Running Variant Quality Score Recalibration for SNPs..."
        local cmd_recal="gatk VariantRecalibrator \
            -R ${REFERENCE_GENOME} \
            -V ${vcf_input} \
            --resource:hapmap,known=false,training=true,truth=true,prior=15.0 ${HAPMAP} \
            --resource:omni,known=false,training=true,truth=false,prior=12.0 ${OMNI} \
            --resource:1000G,known=false,training=true,truth=false,prior=10.0 ${ONEKGENOME} \
            --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ${DBSNP} \
            -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
            -mode SNP \
            --max-gaussians 4 \
            -tranche 100.0 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 95.0 -tranche 90.0 \
            -O ${recal_file} \
            --tranches-file ${tranches_file}"

        #--rscript-file ${rscript_file}"

        local cmd_apply="gatk ApplyVQSR \
            -R ${REFERENCE_GENOME} \
            -V ${vcf_input} \
            --ts-filter-level ${filt_thresh} \
            --tranches-file ${tranches_file} \
            --recal-file ${recal_file} \
            -mode SNP \
            -O ${recal_vcf}"

        #-an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
        execute_and_time "VQSR Recalibration" "Running VQSR recalibration for SNPs." "$cmd_recal"
        execute_and_time "VQSR Application" "Applying VQSR for SNPs." "$cmd_apply"
        #echo "Recalibrated VCF generated at: ${recal_vcf}"
    fi
}


filter_variants() {
    local reference_genome="$1"
    local input_vcf="$2"
    local output_dir=$(dirname "$input_vcf") 
    local outdir="$output_dir/filt_results"
    mkdir -p $outdir
    local output_vcf=$outdir/"$(basename $input_vcf .vcf).varfilter.vcf"
    local filter_name="Low_depth10"
    local filter_expression="DP < 10"
    echo "$output_vcf"
    if [[ ! -s "${output_vcf}" ]]; then
        local cmd="gatk VariantFiltration \
            -R $reference_genome \
            -V $input_vcf \
            -O $output_vcf \
            --filter-name \"$filter_name\" \
            --filter-expression \"$filter_expression\""

        execute_and_time "Filter Variants" "Applying filters to variants based on quality criteria." "$cmd"
        #echo "Filtered VCF created at: $output_vcf"
    fi
}


generate_pass_variants_vcf() {
    local input_vcf="$1"
    local output_dir=$(dirname "$input_vcf") 
    local outdir="$output_dir"
    local output_vcf=$outdir/"$(basename $input_vcf .vcf).pass.vcf"
    echo "$output_vcf"
    if [[ ! -s "${output_vcf}" ]]; then
        # Filter for PASS variants and create a VCF
        local cmd_filter="bcftools view -f 'PASS,.' -O v -o $output_vcf $input_vcf"
        execute_and_time "Extract PASS Variants" "Extracting PASS variants to create analysis-ready VCF." "$cmd_filter"

        # Compress the new VCF file
        local cmd_compress="bgzip -c $output_vcf > ${output_vcf}.gz"
        execute_and_time "Compress VCF" "Compressing the VCF file." "$cmd_compress"

        # Index the compressed VCF file
        local cmd_index="tabix -p vcf ${output_vcf}.gz"
        execute_and_time "Index VCF" "Indexing the compressed VCF file." "$cmd_index"
        #echo "Analysis-ready VCF created and indexed at: ${output_vcf}.gz"
    fi
}


# Function to normalize VCF files using bcftools to remove duplicate records
normalize_dedup_vcf() {
    local input_vcf="$1"
    local output_dir=$(dirname "$input_vcf") 
    local outdir="$output_dir"
    local output_vcf=$outdir/"$(basename $input_vcf .vcf).norm.vcf"
    echo "$output_vcf"

    # Define the normalization command
    local cmd="bcftools norm -D -o ${output_vcf} -O v ${input_vcf}"

    # Execute the command and time it
    execute_and_time "Normalize VCF" "Normalizing VCF file and removing duplicate records to retain the highest quality variants." "$cmd"

    # Output the path to the normalized VCF
    #echo "Normalized VCF generated at: ${full_output_vcf}"
}


# Function to run BEAGLE for all chromosomes
immpute_all_chromosomes() {
    local vcf_file="$1"
    local ref_panel="$2"
    local genetic_map_dir="$3"
    local output_dir="$4"
    local chr_prefix="$5"  # Boolean to include 'chr' prefix in chromosome names
    local outdir="${output_dir}/impute_results"
    mkdir -p "${outdir}"

    # Phase and impute all chromosomes (1-22 and X)
    for chromosome in {1..22} X; do
        local chr_name=$chromosome
        if [ "$chr_prefix" = true ]; then
            chr_name="chr${chromosome}"
        fi

        local map_file="${genetic_map_dir}/plink.${chr_name}.GRCh37.map"
        if [ ! -f "$map_file" ]; then
            echo "Genetic map file not found: $map_file"
            continue  # Skip to next chromosome if map file is missing
        fi

        local out_path="${outdir}/${vcf_file##*/}_chr${chr_name}"

        # Construct the BEAGLE command
        local cmd="java -Xmx1024m -jar beagle.22Jul22.46e.jar \
            chrom=$chr_name \
            gt=$vcf_file \
            ref=$ref_panel \
            map=$map_file \
            out=$out_path \
            impute=true \
            nthreads=10 \
            gp=true"

        # Execute the command with timing
        execute_and_time "BEAGLE Imputation Chromosome $chr_name" "Imputing genetic data for chromosome $chr_name using BEAGLE." "$cmd"
    done
}


# Function to run BEAGLE for all chromosomes with imputation control
phase_impute_all_chromosomes() {
    local vcf_file="$1"
    local ref_panel="$2"
    local genetic_map_dir="$3"
    local output_dir="$4"
    local chr_prefix="$5"  # Boolean to include 'chr' prefix in chromosome names
    local impute="$6"     # Boolean to control imputation

    local result_folder="impute_results"
    if [ "$impute" = false ]; then
        result_folder="phased_results"
    fi

    local outdir="${output_dir}/${result_folder}"
    mkdir -p "${outdir}"

    # Phase and/or impute all chromosomes (1-22 and X)
    for chromosome in {1..22} X; do
        local chr_name=$chromosome
        if [ "$chr_prefix" = true ]; then
            chr_name="chr${chromosome}"
        fi

        local map_file="${genetic_map_dir}/plink.${chr_name}.GRCh37.map"
        if [ ! -f "$map_file" ]; then
            echo "Genetic map file not found: $map_file"
            continue  # Skip to next chromosome if map file is missing
        fi

        local file_suffix="_chr${chr_name}"
        if [ "$impute" = false ]; then
            file_suffix="_phased_chr${chr_name}"
        fi
        local out_path="${outdir}/${vcf_file##*/}${file_suffix}"

        # Construct the BEAGLE command
        local cmd="java -Xmx1024m -jar beagle.22Jul22.46e.jar \
            chrom=$chr_name \
            gt=$vcf_file \
            ref=$ref_panel \
            map=$map_file \
            out=$out_path \
            impute=$impute \
            nthreads=10 \
            gp=true"

        # Execute the command with timing
        local action="Imputing"
        if [ "$impute" = false ]; then
            action="Phasing"
        fi
        execute_and_time "BEAGLE $action Chromosome $chr_name" "$action genetic data for chromosome $chr_name using BEAGLE." "$cmd"
    done
}


# Function to run BEAGLE for all chromosomes with expanded imputation control
boolean_phase_impute_all_chromosomes() {
    local vcf_file="$1"
    local ref_panel="$2"
    local genetic_map_dir="$3"
    local output_dir="$4"
    local chr_prefix="$5"  # Boolean to include 'chr' prefix in chromosome names
    local impute_option="$6" # Can be true, false, or both

    local impute_actions=("true")
    if [ "$impute_option" = "false" ]; then
        impute_actions=("false")
    elif [ "$impute_option" = "both" ]; then
        impute_actions=("true" "false")
    fi

    for impute in "${impute_actions[@]}"; do
        local result_folder="impute_results"
        if [ "$impute" = "false" ]; then
            result_folder="phased_results"
        fi
        local outdir="${output_dir}/${result_folder}"
        mkdir -p "${outdir}"

        # Phase and/or impute all chromosomes (1-22 and X)
        for chromosome in {1..22} X; do
            local chr_name=$chromosome
            if [ "$chr_prefix" = true ]; then
                chr_name="chr${chromosome}"
            fi

            local map_file="${genetic_map_dir}/plink.${chr_name}.GRCh37.map"
            if [ ! -f "$map_file" ]; then
                echo "Genetic map file not found: $map_file"
                continue  # Skip to next chromosome if map file is missing
            fi

            local file_suffix="_chr${chr_name}"
            if [ "$impute" = "false" ]; then
                file_suffix="_phased_chr${chr_name}"
            fi
            local out_path="${outdir}/${vcf_file##*/}${file_suffix}"

            # Construct the BEAGLE command
            if [[ ! -s "$output_vcf" ]]; then
                local cmd="beagle \
                    chrom=$chr_name \
                    gt=$vcf_file \
                    ref=$ref_panel \
                    map=$map_file \
                    out=$out_path \
                    impute=$impute \
                    nthreads=10 \
                    gp=true"

                # Execute the command with timing
                local action="Imputing"
                if [ "$impute" = "false" ]; then
                    action="Phasing"
                fi
                execute_and_time "BEAGLE $action Chromosome $chr_name" "$action genetic data for chromosome $chr_name using BEAGLE." "$cmd"
            fi
        done
    done
}


# BEAGLE Function Description
#   This function executes the process of running the BEAGLE tool for
#   imputation and phasing across all human chromosomes (1 to 22 and X).
#   Provides flexibility in operation to choose between phasing only,
#   imputation only, or both.

# Function to run BEAGLE for all chromosomes with expanded imputation control
run_beagle_for_all_chromosomes() {
    local vcf_file="$1"
    local ref_panel="$2"
    local genetic_map_dir="$3"
    local output_dir="$4"
    local chr_prefix="$5" # Boolean to include 'chr' prefix in chromosome names
    local impute_option="$6" # Can be phase, impute, both (equals false, true and false/true) 
    #local outdir="${output_dir}/$sample/results"
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
        local outdir="${output_dir}/$sample/results/${result_folder}"
        mkdir -p "${outdir}"

        # Process all chromosomes (1-22 and X)
        for chromosome in {1..22} X; do
            local chr_name=$chromosome
            [ "$chr_prefix" = true ] && chr_name="chr${chromosome}" && map_file="${genetic_map_dir}/plink.${chr_name}.GRCh37.map"
            #[ "$chr_prefix" = false ] && local map_file="${genetic_map_dir}/plink.chr${chr_name}.GRCh37.map"
            
            [ "$chr_prefix" = false ] && local map_file="${genetic_map_dir}/plink.chr${chr_name}.GRCh37.map"
            #local map_file="${genetic_map_dir}/plink.chr${chr_name}.GRCh37.map
            echo -e "\nProcessing chromosome: $chr_name"

            if [ ! -f "$map_file" ]; then
                echo "Genetic map file not found: $map_file"
                continue  # Skip missing map files
            fi

            local file_suffix=".chr${chr_name}"
            if [ "$impute" = "false" ]; then
                file_suffix=".phased.chr${chr_name}"
            elif [ "$impute" = "true" ]; then
                file_suffix=".imputed.chr${chr_name}"
            fi

            #local file_suffix="_chr${chr_name}${impute:0:1}"
            local vcf_outbase="${outdir}/${vcf_file##*/}${file_suffix}"
            local vcf_outbase="${outdir}/$(basename ${vcf_file} .vcf)${file_suffix}"
            local output_vcf=${vcf_outbase}.vcf.gz
            echo -e "\nOUTPUT BEAGLE VCF: $output_vcf"
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


# Function to ensure training data has necessary annotations
ensure_annotated_training_data() {
    local training_vcf="$1"
    local annotated_vcf="$2"
    echo "Ensuring training data ${training_vcf} has necessary annotations..."
    gatk VariantAnnotator \
        -R ${REFERENCE_GENOME} \
        -V ${training_vcf} \
        -O ${annotated_vcf} \
        -A QualByDepth \
        -A MappingQuality \
        -A MappingQualityRankSumTest \
        -A ReadPosRankSumTest \
        -A FisherStrand \
        -A StrandOddsRatio
}


# Function for Variant Quality Score Recalibration (VQSR) for SNPs
vqsr_annotated_version() {
    local sample="$1"
    local vcf_input="$2"
    local output_dir="$3"
    local outdir="${output_dir}/${sample}/results"
    local recal_vcf="${outdir}/${sample}_VQSR.snp.recal.vcf"

    mkdir -p "${outdir}"  # Ensure the output directory exists

    # Ensure annotated versions of the training data are available
    ensure_annotated_training_data "${HAPMAP}" "${outdir}/hapmap_annotated.vcf"
    ensure_annotated_training_data "${OMNI}" "${outdir}/omni_annotated.vcf"
    ensure_annotated_training_data "${ONEKGENOME}" "${outdir}/1000G_annotated.vcf"
    ensure_annotated_training_data "${DBSNP}" "${outdir}/dbsnp_annotated.vcf"

    # VQSR command setup
    local cmd_recal="gatk VariantRecalibrator \
        -R ${REFERENCE_GENOME} \
        -V ${vcf_input} \
        --resource:hapmap,known=false,training=true,truth=true,prior=15.0 ${outdir}/hapmap_annotated.vcf \
        --resource:omni,known=false,training=true,truth=false,prior=12.0 ${outdir}/omni_annotated.vcf \
        --resource:1000G,known=false,training=true,truth=false,prior=10.0 ${outdir}/1000G_annotated.vcf \
        --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ${outdir}/dbsnp_annotated.vcf \
        -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
        -mode SNP \
        --max-gaussians 4 \
        -tranche 100.0 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 95.0 -tranche 90.0 \
        -O ${recal_file} \
        --tranches-file ${tranches_file} \
        --rscript-file ${rscript_file}"

    local cmd_apply="gatk ApplyVQSR \
        -R ${REFERENCE_GENOME} \
        -V ${vcf_input} \
        --ts-filter-level 99.5 \
        --tranches-file ${tranches_file} \
        --recal-file ${recal_file} \
        -mode SNP \
        -O ${recal_vcf}"
        
    #--ts-filter-level 99.9 \
    execute_and_time "VQSR Recalibration" "Running VQSR recalibration for SNPs." "$cmd_recal"
    execute_and_time "VQSR Application" "Applying VQSR for SNPs." "$cmd_apply"
    echo "Recalibrated VCF generated at: ${recal_vcf}"
}


# Assuming execute_and_time is defined elsewhere, here is a placeholder definition
#execute_and_time() {
#    local description="$1"
#    local detailed_description="$2"
#    local command="$3"
#    echo "Starting $description: $detailed_description"
#    local start_time=$(date +%s)
#    eval $command
#    local status=$?
#    local end_time=$(date +%s)
#    local elapsed_time=$((end_time - start_time))
#    if [ $status -ne 0 ]; then
#        echo "Error during $description. Exiting with status $status."
#        exit $status
#    else
#        echo "$description completed successfully in $elapsed_time seconds."
#    fi
#}


# Function to perform various BAM file analyses
qcanalyze_bam_file() {
    local sample="$1"
    local bam_file="$2"
    local output_dir="$3"
    local reference_fasta="${REFERENCE_GENOME}"
    
    # Ensure output directory exists
    outdir=$output_dir/$sample/qc_results
    mkdir -p "$outdir"
    
    # Qualimap BAM quality check
    qualimap_file="$outdir/${sample}_qualimap_report.pdf"
    if [[ ! -s "$qualimap_file" ]]; then 
        execute_and_time "Qualimap BAM QC" "Running Qualimap on the BAM file." \
        "qualimap bamqc --java-mem-size=50G -bam $bam_file -outdir $outdir -outfile ${sample}_qualimap_report.pdf -outformat PDF"
    fi 

    # Picard CollectAlignmentSummaryMetrics
    #picard_file="$outdir/${sample}_picard_output_metrics.txt"
    #if [[ ! -s "$picard_file" ]]; then
    #    execute_and_time "Picard CollectAlignmentSummaryMetrics" "Collecting alignment summary metrics." \
    #    "picard CollectAlignmentSummaryMetrics R=$reference_fasta I=$bam_file O=${outdir}/${sample}_picard_output_metrics.txt"
    #fi
    
    # Picard CollectInsertSizeMetrics
    picard_file2="${outdir}/${sample}_picard_insert_size_histogram.pdf"
    if [[ ! -s "$picard_file2" ]]; then
        execute_and_time "Picard CollectInsertSizeMetrics" "Collecting insert size metrics." \
        "picard -Xmx50g CollectInsertSizeMetrics I=$bam_file O=${outdir}/picard_insert_metrics.txt H=${outdir}/${sample}_picard_insert_size_histogram.pdf"
    fi
    
    # Samtools statistics
    samtool_file="${outdir}/${sample}_samtools.fullstats"
    if [[ ! -s "$samtool_file" ]]; then
        execute_and_time "Samtools Stats" "Generating statistics with Samtools." \
        "samtools stats $bam_file > ${outdir}/${sample}_samtools.fullstats"
    fi
    
    # Samtools statistics
    samtool_file2="${outdir}/${sample}_samtools.flagstats"
    if [[ ! -s "$samtool_file2" ]]; then
        execute_and_time "Samtools Stats" "Generating statistics with Samtools." \
        "samtools flagstat $bam_file > ${outdir}/${sample}_samtools.flagstats"
    fi
    
    # Bedtools genome coverage
    bedtool_file="${outdir}/${sample}_bedtool_coverage.txt"
    if [[ ! -s "$bedtool_file" ]]; then
        execute_and_time "Bedtools Genome Coverage" "Calculating genome coverage with Bedtools." \
        "bedtools genomecov -ibam $bam_file -dz > ${outdir}/${sample}_bedtool_coverage.txt"
    fi
}


# Function to generate QC reports for sequencing data
generate_qc_reports() {
    local sample="$1"
    local fastq1="$2"
    local fastq2="$3"
    local output_dir="$4"  # Ensure output_dir is passed as a parameter
    local qc_output_dir="${output_dir}/QC_reports/${sample}"

    # Function to create a directory if it does not exist
    create_dir_if_not_exists() {
        local dir="$1"
        if [[ ! -d "$dir" ]]; then
            mkdir -p "$dir"
        fi
    }

    # Create QC output directory if it does not exist
    create_dir_if_not_exists "$qc_output_dir"

    # Command to generate QC reports using FastQC
    cmd="fastqc -o ${qc_output_dir} ${fastq1} ${fastq2}"
    execute_and_time "FastQC Analysis" "Generating FastQC reports for raw fastq data." "$cmd"

    # Command to aggregate QC reports using MultiQC
    #local multiqc_cmd="multiqc -o ${qc_output_dir} ${qc_output_dir}"
    #execute_and_time "MultiQC Analysis" "Aggregating QC reports using MultiQC." "$multiqc_cmd"    
}


# Main script logic
#if [[ -z "$genome_fasta" || -z "$chromosome" || -z "$fastq1" || -z "${genome-index}" ]]; then
#    echo "Error: Missing required arguments."
#    usage
#    exit 1
#fi

# Main script logic
if [[ "$genome_index" != "true" ]]; then
    # Check if essential arguments are missing when not in genome-index mode
    if [[ -z "$genome_fasta" || -z "$chromosome" || -z "$fastq1" ]]; then
        echo "Error: Missing required arguments."
        usage
        exit 1
    fi
else
    # Check for genome_fasta and chromosome when in genome-index mode
    if [[ -z "$genome_fasta" || -z "$chromosome" ]]; then
        echo "Error: Missing required arguments --genome-fasta or --chromosome."
        usage
        exit 1
    fi
fi


# Check tools availability
if ! command -v samtools &> /dev/null || ! command -v bowtie2 &> /dev/null || ! command -v bwa &> /dev/null; then
    echo "Error: Required tools (samtools, bowtie2, bwa) are not available."
    exit 1
fi


main_function(){
    local version=1.0.0
    echo -e "\n\n********PRIMARY SECTION OF SC-VARIANT CALLING PIPELINE: V$version********\n\n"
    FILE_DIR=${output_dir}
    output_dir="${output_dir}/genome_indexes"
    create_dir_if_not_exists "${output_dir}"
    
    # Handle genome indexing and optional chromosome extraction based on --genome-index flag
    if [[ "$genome_index" == "true" ]]; then
        echo "Starting genome indexing only..."
        genomefasta=$(index_genome "$genome_fasta" "$output_dir")
        echo -e "Indexed genome fasta file: $genomefasta\n\n"

        # Check if the chromosome flag is set to 'all'
        if [[ "$chromosome" == "all" ]]; then
            echo "Processing entire genome..."
            ref_index=$(build_reference_index "$genomefasta" "$aligner")
        else
            chr_list=$(echo "$chromosome" | tr ',' ' ')
            echo "Extracting chromosome to create new genomefasta: ${chr_list}"
            genomefasta=$(extract_chromosome "$genomefasta" "$chr_list" "$output_dir")
            echo -e "\n\nChromosome fasta filename: $genomefasta \n\n"

            echo "Building reference index for chromosome..."
            ref_index=$(build_reference_index "$genomefasta" "$aligner")
        fi

        echo -e "\n\n***********GENOME-INDEXING/EXTRACTION COMPLETED***********"
        echo -e "\n\nReference index filename: $ref_index"
        exit 0 # exit successfully

    else
        # Execute workflow
        echo "Starting genome indexing..."
        genomefasta=$(index_genome "$genome_fasta" "$output_dir")
        echo -e "Indexed genome fasta file: $genomefasta\n\n"

        # Check if the chromosome flag is set to 'all'
        if [[ "$chromosome" == "all" ]]; then
            echo "Processing entire genome..."
            ref_index=$(build_reference_index "$genomefasta" "$aligner")
        else
            chr_list=$(echo "$chromosome" | tr ',' ' ')
            echo "Extracting chromosome to create new genomefasta: ${chr_list}"
            genomefasta=$(extract_chromosome "$genomefasta" "$chr_list" "$output_dir")
            echo -e "\n\nChromosome fasta filename: $genomefasta \n\n"

            echo "Building reference index for chromosome..."
            ref_index=$(build_reference_index "$genomefasta" "$aligner")
        fi
        echo -e "\n\nReference index filename: $ref_index"

        # To be consistent with some functions variable names
        REFERENCE_GENOME=${genomefasta}
        output_dir="${FILE_DIR}/sample_outputs"
        create_dir_if_not_exists "${output_dir}"
       
        #echo -e "\nProcessing fastqc report for sample: $sample"
        #generate_qc_reports "$sample" "$fastq1" "$fastq2" "$output_dir"
        #echo -e "\nFastqc report for sample: $sample completed ..."

        # Read each line in the input file and apply parse_fastq_headers
        parsed_rgid=$(parse_fastq_headers $fastq1)
        echo -e "\nParsed Uniqe RG ID: $parsed_rgid"

        bam_file=$(align_reads "$ref_index" "$fastq1" "$fastq2" "$aligner" "$output_dir" "$parsed_rgid")
        echo -e "\n\nBam filename: $bam_file"

        sorted_bam=$(sort_bam "$sample" "$bam_file" "$output_dir")
        echo -e "\n\nSorted bam filename: $sorted_bam"

        dedup_bam=$(mark_duplicates "$sample" "$sorted_bam" "$output_dir")
        echo -e "\n\nDedup bam filename: $dedup_bam"

        # Perform post-base quality recalibration
        #post_base_quality_recalibration "$sample" "$recal_bam" "$output_dir"
        #echo -e "\n\nPost base quality score recalibration and visualization: $recal_bam"

        #qcanalyze_bam_file=$(qcanalyze_bam_file "${sample}" "$dedup_bam" "$output_dir")
        #echo -e "\n\nQC deduped bam_file: $dedup_bam"

        create_fasta_dictionary "$genomefasta"
        echo -e "Created SequenceDictionary: $genomefasta\n\n"

        recompress_and_index_vcf "${DBSNP}"
        echo -e "Indexed vcf file: ${DBSNP}\n\n"

        recompress_and_index_vcf "${HAPMAP}"
        echo -e "Indexed vcf file: ${HAPMAP}\n\n"

        recompress_and_index_vcf "${OMNI}"
        echo -e "Indexed vcf file: ${OMNI}\n\n"

        recompress_and_index_vcf "${ONEKGENOME}"
        echo -e "Indexed vcf file: ${ONEKGENOME}\n\n"

        ## genome fasta file is equivalent to REFERENCE_GENOME file
        ##check_and_conditionally_reheader_vcf "$genome_fasta" "$DBSNP" "$CONTIGS_CONVERSION_FILE"
        ##echo -e "Check and conditionally reheader vcf file: ${DBSNP}\n\n"
        ##check_and_conditionally_rename_vcf "$genome_fasta" "$DBSNP" "$CONTIGS_CONVERSION_FILE"
        ##echo -e "Check and conditionally rename vcf file: ${DBSNP}\n\n"

        recal_bam_full=$(base_quality_recalibration "$sample" "$dedup_bam" )
        echo -e "\n\nRecalibrated bam_file: $recal_bam"

        gvcfs=$(variant_calling "$sample" "$recal_bam" )
        echo -e "\n\nGenotype calling from GVCFs: $gvcfs"

        genotype_vcf=$(genotype_gvcf "$sample" "$gvcfs"  "single")
        ##genotype_vcf=$(genotype_gvcf "$sample" "$gvcfs" "$output_dir" "multi")
        echo -e "\n\nGenotype final output file: $genotype_vcf"

        ##vqsr_recal_vcf=$(variant_quality_score_recalibration "$sample" "$genotype_vcf" "$output_dir")
        ##echo -e "\n\nVQSR output file: $vqsr_recal_vcf"
        vqsr_recal_vcf_full=$(variant_quality_score_recalibration "$sample" "$genotype_vcf")
        echo -e "\n\nVQSR output file: $vqsr_recal_vcf"

        filt_variant_vcf=$(filter_variants "$genomefasta" "$vqsr_recal_vcf")
        echo -e "\n\nVQSR variant vcf file filtered by depth: $filt_variant_vcf"

        passed_variant_vcf=$(generate_pass_variants_vcf "$filt_variant_vcf")
        echo -e "\n\nVQSR variant vcf file with pass status: $passed_variant_vcf"

        ref_panel="$REFERENCE_PANEL"
        normalized_refpanel_vcf=$(normalize_dedup_vcf "$ref_panel")
        echo "Normalized refpanel deduped vcf generated at: ${normalized_refpanel_vcf}"
      
        ##beagle_outdir=$(run_beagle_for_all_chromosomes "$vcf_file" "$ref_panel" "$genetic_map_dir" "$output_dir" false false)
        #beagle_outdir=$(run_beagle_for_all_chromosomes "$normalized_dedup_vcf" "$ref_panel" "$genetic_map_dir" "$output_dir" false "both")
        genetic_map_dir="$GENETIC_MAP"
        beagle_outdir=$(run_beagle_for_all_chromosomes "$passed_variant_vcf" "$normalized_refpanel_vcf" "$genetic_map_dir" "$output_dir" false "both")
        echo -e "\n\nBeagle ouput vcf files in dir: $beagle_outdir"

        # QC evaluations        
        echo -e "\nProcessing fastqc report for sample: $sample"
        generate_qc_reports "$sample" "$fastq1" "$fastq2" "$output_dir"
        echo -e "\nFastqc report for sample: $sample completed ..."

        # Perform post-base quality recalibration
        post_base_quality_recalibration "$sample" "$recal_bam" "$output_dir"
        echo -e "\n\nPost base quality score recalibration and visualization: $recal_bam"
        
        qcanalyze_bam_file=$(qcanalyze_bam_file "${sample}" "$dedup_bam" "$output_dir")
        echo -e "\n\nQC deduped bam_file: $dedup_bam"
        
   fi         
}


# Call genetic variants and generate QC reports
version=1.0.0
echo -e "\n\n********SC-VARIANT CALLING PIPELINE: V$version********\n\n"
echo -e "\n\nParameters:
            genomefasta : $genome_fasta
            fastq1: $fastq1
            fastq2: $fastq2
            chrom: $chromosome
            sample: $sample
            aligner: $aligner
            output dir: $output_dir
        "

# Execute main_function
echo -e "\nPipeline run starting for sample: $sample"
main_function
echo -e "\n\n******* Process complete. Outputs are in $output_dir *******\n\n"


