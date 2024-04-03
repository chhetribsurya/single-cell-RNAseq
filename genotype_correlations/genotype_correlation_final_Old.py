#!/usr/bin/env python

# Python script for final correlation computation
import pandas as pd
import os
import re
from os.path import join
from scipy.stats import pearsonr, spearmanr
import matplotlib.pyplot as plt
import sys

# ref_file = join("{}".format(output_dir), "{}_ref_filtered.csv".format(base_output_file))
# imputed_file = join("{}".format(output_dir), "{}_gencove_harmonized_subsetted.csv".format(base_output_file))

# File paths as command-line arguments
ref_file = sys.argv[1]
imputed_file = sys.argv[2]
outfile_prefix = sys.argv[3]

# Function to convert genotype to numeric value
def genotype_to_numeric(genotype):
    if genotype == '0/0':
        return 0
    elif genotype in ['0/1', '1/0']:
        return 1
    elif genotype == '1/1':
        return 2
    else:
        return None  # Handle unexpected values

# Extract and save the correlations to output file
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
                #print(f'Processing file: {filename}')
                #print(f'Pearson Correlation: {pearson_corr}')
                #print(f'Spearman Correlation: {spearman_corr}\n')
                out_file.write(f'{filename}\t{pearson_corr}\t{spearman_corr}\n')

# Assuming bi-allelic variants : convert reference genotype calls of 0/0, 0/1, 1/1 etc 
# into a comparable numeric value. For example:
# 0/0 -> 0;
# 0/1 and 1/0 -> 1;
# 1/1 -> 2

print("Reading reference file ...")
ref_genotypes = pd.read_csv(ref_file, sep="\t", header=None, names=["chromosome", "position", "GT"])

print("Reading imputed file ...")
imputed_genotypes = pd.read_csv(imputed_file, sep="\t", header=None, names=["chromosome", "position", "GT"])

print("Merging of reference and imputed file ...")
merged_df = pd.merge(ref_genotypes, imputed_genotypes, on=["chromosome", "position"])

# Rename columns
merged_df.rename(columns={'GT_x': 'GT_ref', 'GT_y': 'GT_imputed'}, inplace=True)

# Apply the conversion
merged_df['GT_ref_numeric'] = merged_df['GT_ref'].apply(genotype_to_numeric)
merged_df['GT_imputed_numeric'] = merged_df['GT_imputed'].apply(genotype_to_numeric)

# Compute Pearson and Spearman correlations
print("Computing correlation ...")
pearson_corr, _ = pearsonr(merged_df['GT_ref_numeric'], merged_df['GT_imputed_numeric'])
spearman_corr, _ = spearmanr(merged_df['GT_ref_numeric'], merged_df['GT_imputed_numeric'])

print(f"\nPearson Correlation: {pearson_corr:.2f}")
print(f"Spearman Correlation: {spearman_corr:.2f}\n")

# Plotting
print("Plotting correlation plot ...")
plt.figure(figsize=(10, 6))
plt.scatter(merged_df['GT_ref_numeric'], merged_df['GT_imputed_numeric'], alpha=0.6)
plt.xlabel('GT_ref_numeric')
plt.ylabel('GT_imputed_numeric')
plt.title(f'Genotype Comparison\nPearson Correlation: {pearson_corr:.2f}, Spearman Correlation: {spearman_corr:.2f}')
plt.grid(True)

# Define the directory for plot files
plot_files_dir = './plot_files'

# Create the plot_files directory if it doesn't exist
if not os.path.exists(plot_files_dir):
    os.makedirs(plot_files_dir)

# Save the plot in the plot_files directory
plot_file_path = os.path.join(plot_files_dir, outfile_prefix + '_genotype_comparison_plot.pdf')
plt.savefig(plot_file_path)

# Save the plot
# plt.savefig(outfile_prefix + '_genotype_comparison_plot.pdf')

# Define the directory and output file
log_files_dir = './log_files'
output_tsv = './correlation_results.tsv'  # Replace with output path

# Process the log files and output correlation results
process_log_files(log_files_dir, output_tsv)
