#!/bin/bash

set -euox pipefail

# Function to run BLAST
run_blast() {
    local query=$1
    local fasta=$2

    local ourbasename=$(basename "$fasta" .fasta)

    # Create BLAST database
    makeblastdb -in $fasta -dbtype nucl -out "$ourbasename.db"

    # Run BLAST
    blastn -max_target_seqs 100000000 -query "$query" -db "$ourbasename.db" -task "blastn-short" -outfmt 6 -evalue 1 -word_size 7 -perc_identity 90 > "$ourbasename.blast.tsv"

    echo "BLAST output generated: $ourbasename.blast.tsv"
}


convert_fastq_to_fasta() {
    local path="$1"
    local base_path=$(basename "$path" .fastq.gz)
    local fasta_path="${base_path}.fasta"
    gzcat "$path" | seqtk seq -a > "$fasta_path"
    echo "$fasta_path"
}


# Light chain

fasta_path=$(convert_fastq_to_fasta ../input/outname_AV_IgYKSTOP_8A_S2_R1_001_atleast-2.fastq.gz)
run_blast around-deleted-8a-11a.fasta $fasta_path

fasta_path=$(convert_fastq_to_fasta ../input/outname_AV_IgYKSTOP_9A_S3_R1_001_atleast-2.fastq.gz)
run_blast around-deleted-9a.fasta $fasta_path

fasta_path=$(convert_fastq_to_fasta ../input/outname_AV_IgYKSTOP_11A_S1_R1_001_atleast-2.fastq.gz)
run_blast around-deleted-8a-11a.fasta $fasta_path


# Heavy chain

for path in ../input/outname_[567]*atleast-2.fastq.gz; do
    fasta_path=$(convert_fastq_to_fasta "$path")
    run_blast around-inserted-igh.fasta $fasta_path
done
