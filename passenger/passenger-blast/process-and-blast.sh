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
    blastn -max_target_seqs 100000000 -query "$query" -db "$ourbasename.db" -task "blastn-short" -outfmt 6 -evalue 1000 -word_size 7 -perc_identity 95 > "$ourbasename.blast.tsv"

    echo "BLAST output generated: $ourbasename.blast.tsv"
}


# Light chain

bash ./filter-to-atleast2.sh ../input/AV_IgYKSTOP_8A_collapse-unique.fastq.gz
run_blast around-deleted-8a.fasta AV_IgYKSTOP_8A_collapse-unique.atleast2.fasta

bash ./filter-to-atleast2.sh ../input/AV_IgYKSTOP_9A_collapse-unique.fastq.gz
run_blast around-deleted-9a.fasta AV_IgYKSTOP_9A_collapse-unique.atleast2.fasta

bash ./filter-to-atleast2.sh ../input/AV_IgYKSTOP_11A_collapse-unique.fastq.gz
run_blast around-deleted-11a.fasta AV_IgYKSTOP_11A_collapse-unique.atleast2.fasta


# Heavy chain

for path in ../input/outname_[567]*atleast-2.fastq.gz; do
    base_path=$(basename "$path" .fastq.gz)
    fasta_path="$base_path".fasta
    gzcat "$path" | seqtk seq -a > $fasta_path
    run_blast around-inserted-igh.fasta $fasta_path
done
