import os
import gzip

def count_reads_fastq_gz(file_path):
    count = 0
    with gzip.open(file_path, 'rt') as f:  # Open the gzipped fastq file
        for line in f:
            if line.startswith('@'):  # Each read starts with '@'
                count += 1
    return count

# Directory containing the fastq.gz files
directory = '/data/abattle4/lakshmi/cuomo_2020/pipeline_files/concat_files'

# Output file to save the read counts
output_file = './readcount_outputfile.txt'

# Print and count only the R1 fastq.gz files in the directory
print("Processing Fastq.gz files in the directory:")
with open(output_file, 'w') as out:
    for filename in os.listdir(directory):
        if filename.endswith('_R1.fastq.gz'):
            print(f"\nProcessing file: {filename} ...")
            file_path = os.path.join(directory, filename)
            read_count = count_reads_fastq_gz(file_path)
            print(f"{filename}: {read_count} reads")
            out.write(f"{filename}: {read_count} reads\n")

