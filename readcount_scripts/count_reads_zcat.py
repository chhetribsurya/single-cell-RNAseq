import os
import subprocess

def count_reads_unix(file_path):
    command = f"zcat {file_path} | wc -l"
    process = subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

    if process.returncode == 0:
        line_count = int(process.stdout.strip())
        return line_count // 4
    else:
        print(f"Error processing file {file_path}: {process.stderr}")
        return None

# Directory containing the fastq.gz files
directory = '/data/abattle4/lakshmi/cuomo_2020/pipeline_files/concat_files'

# Output file to save the read counts
output_file = './readcount_zcatoutput.txt'

# Print and count only the R1 fastq.gz files in the directory
print("Processing Fastq.gz files in the directory:")
with open(output_file, 'w') as out:
    for filename in os.listdir(directory):
        if filename.endswith('_R1.fastq.gz'):
            print(f"\nProcessing file: {filename} ...")
            file_path = os.path.join(directory, filename)
            read_count = count_reads_unix(file_path)
            if read_count is not None:
                print(f"{filename}: {read_count} reads")
                out.write(f"{filename}: {read_count} reads\n")

