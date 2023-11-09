# for file in ../../nextflow/data/passenger/test_run/unique_with_at_least_two_reads/*.fastq
# for file in ../../nextflow/data/passenger/8A/AV_IgYKSTOP_8A_UniqueSequences_atleast-2.fastq
# query="around-deleted-corrected.fasta"
query="around-deleted-9a.fasta"
for file in ../../nextflow/data/passenger/9A/9A_collapse-unique_atleast-2.fastq
do
    filename=$(basename "$file" .fastq)
    seqtk seq -A "$file" > "$filename.fasta"
    makeblastdb -in "$filename.fasta" -dbtype nucl -out "$filename.db"
    blastn -max_target_seqs 100000000 -query "$query" -db "$filename.db" -task "blastn-short" -outfmt 6 -evalue 1000 -word_size 7 -perc_identity 95 > "$filename.blast.tsv"
done
