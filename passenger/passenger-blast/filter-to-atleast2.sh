#!/bin/bash

# Exit when any command fails and catch errors in pipelines
set -euo pipefail

# Check if the file path is provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <path_to_fastq_gz>"
    exit 1
fi

# Input gzipped FASTQ file
input_fastq_gz="$1"

# Check if the file exists
if [ ! -f "$input_fastq_gz" ]; then
    echo "File not found: $input_fastq_gz"
    exit 1
fi

# Construct output FASTA filename
output_fasta=$(basename "$input_fastq_gz" | sed 's/\.fastq\.gz$/.atleast2.fasta/')

# Create a temporary file for sequence IDs
temp_seq_ids=$(mktemp)

# Extract sequence IDs that do not contain 'CONSCOUNT=1'
gunzip -c "$input_fastq_gz" | grep -v 'CONSCOUNT=1' | grep '^@' | sed 's/^@//' > "$temp_seq_ids"

# Use seqtk to subset the original FASTQ file and convert to FASTA
seqtk subseq "$input_fastq_gz" "$temp_seq_ids" | seqtk seq -a > "$output_fasta"

# Remove temporary file
rm "$temp_seq_ids"

echo "Output generated: $output_fasta"
