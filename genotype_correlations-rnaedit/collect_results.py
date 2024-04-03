import os
import re

def extract_correlations(file_path):
    """Extracts Pearson and Spearman correlations from a given log file."""
    pearson_corr = spearman_corr = None
    with open(file_path, 'r') as file:
        for line in file:
            if 'Pearson Correlation:' in line:
                pearson_corr = re.search(r'Pearson Correlation: ([0-9.]+)', line).group(1)
            elif 'Spearman Correlation:' in line:
                spearman_corr = re.search(r'Spearman Correlation: ([0-9.]+)', line).group(1)
    return pearson_corr, spearman_corr


def process_log_files(directory, output_file):
    """Processes all log files in the given directory and writes results to a TSV file."""
    with open(output_file, 'w') as out_file:
        out_file.write('Filename\tPearson Correlation\tSpearman Correlation\n')
        for filename in os.listdir(directory):
            if filename.endswith('.log'):  # Assuming log files end with '.log'
                file_path = os.path.join(directory, filename)
                pearson_corr, spearman_corr = extract_correlations(file_path)

                # Print the processing file and correlations
                print(f'Processing file: {filename}')
                print(f'Pearson Correlation: {pearson_corr}')
                print(f'Spearman Correlation: {spearman_corr}\n')

                out_file.write(f'{filename}\t{pearson_corr}\t{spearman_corr}\n')


# Define the directory and output file
#log_files_dir = '/data/abattle4/surya/datasets/for_lakshmi/final_scripts/log_files'
log_files_dir = './log_files'
output_tsv = './correlation_results.tsv'  # Replace with your desired output path

# Process the log files
process_log_files(log_files_dir, output_tsv)

