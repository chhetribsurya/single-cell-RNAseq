#!/bin/bash

#_____________________________________________________CONFIG FILE__________________________________________________________#

# --- CSV File Setup ---
CSV_FILE="sample_fastq_paths.csv"
CSV_FILE="sample_fastq_paths-test.csv"


# --- Path Configuration ---
CHROMOSOME="all"
CHROMOSOME="1,2,17,20"
CHROMOSOME="1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X"
GENOME_FASTA="/scratch16/abattle4/surya/datasets/for_gatk/resource_bundles/GRCh37/hs37d5.fa"
OUTPUT_DIR="/data/abattle4/surya/datasets/for_lakshmi/variant_calling/variant_call_STAR"
SCRIPT_DIR="/data/abattle4/surya/datasets/for_lakshmi/variant_calling/scripts"


# --- Environment Setup ---
module load r/4.0.2
mkdir -p "$OUTPUT_DIR"


# Define aligner
ALIGNER="bowtie2"
ALIGNER="bwa"
ALIGNER="STAR"


# Perform genome-indexing
bash ${SCRIPT_DIR}/align-final-edit.sh \
    --genome-fasta $GENOME_FASTA \
    --chromosome ${CHROMOSOME} \
    --output-dir $OUTPUT_DIR \
    --aligner STAR \
    --genome-index


# Remove temp dir
if [[ -d "${OUTPUT_DIR}/resultOutfile" ]]; then
    rm -r "${OUTPUT_DIR}/resultOutfile"
fi


echo -e "\n\n********PRE-INDEXING OF GENOME COMPLETED*********"
echo -e "\nNow proceeding with variant calling pipeline for all samples ..."


# --- Create Task File from CSV for sample submission---
TASK_FILE="${OUTPUT_DIR}/script_tasks.txt"
echo "" > "$TASK_FILE"


tail -n +2 "$CSV_FILE" | while IFS=',' read -r sample fastq1 fastq2; do
    # Construct the command for each sample
    echo "bash $SCRIPT_DIR/align-final-edit.sh \
        --genome-fasta $GENOME_FASTA \
        --chromosome $CHROMOSOME \
        --output-dir $OUTPUT_DIR \
        --fastq1 $fastq1 \
        --fastq2 $fastq2 \
        --aligner $ALIGNER \
        --assaytype "RNA-seq" \
        --sample $sample" >> "$TASK_FILE"
done


echo -e "\n\nFind job submitting task file at: $TASK_FILE"

# --- SLURM Batch Script Creation ---
cat <<EOF > ${OUTPUT_DIR}/run_array.sh
#!/bin/bash
#SBATCH --job-name=variant-calling
#SBATCH --output=${OUTPUT_DIR}/slurmlogs/slurm_%A_%a.out
#SBATCH --error=${OUTPUT_DIR}/slurmlogs/slurm_%A_%a.err
#SBATCH --array=2-$(wc -l < $TASK_FILE)%20
#SBATCH --mem=70G
#SBATCH --cpus-per-task=20
#SBATCH --time=72:00:00


# module load
ml r/4.0.2

# conda env
eval "\$(conda shell.bash hook)"
conda activate gatk_env

# Read the command from the task file using SLURM_ARRAY_TASK_ID
CMD=\$(sed -n "\${SLURM_ARRAY_TASK_ID}p" $TASK_FILE)
eval "\$CMD"
EOF


# --- Submit SLURM Job ---
sbatch ${OUTPUT_DIR}/run_array.sh




