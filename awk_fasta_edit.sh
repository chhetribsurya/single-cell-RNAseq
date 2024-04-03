#!/bin/bash

# File path to the input FASTA file
fasta_file="/scratch16/abattle4/surya/worksets/for_lakshmi/final_fasta_edits/hs37d5.fa"

# Output file
output_file="./$(basename $fasta_file .fa).edit.fa"

# Using awk to prepend 'chr' to lines for chromosomes 1 to 22, X, and Y
echo -e "\nProcessing file: $(basename $fasta_file)"

awk '
    /^>1|^>2|^>3|^>4|^>5|^>6|^>7|^>8|^>9|^>10|^>11|^>12|^>13|^>14|^>15|^>16|^>17|^>18|^>19|^>20|^>21|^>22|^>X|^>Y/ {
        sub(/^>/, ">chr");
    }
    { print }
' "$fasta_file" > "$output_file"

# Replace original file with the modified one
#mv "$output_file" "$fasta_file"

echo "FASTA file updated with chr prefix."

