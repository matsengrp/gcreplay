# seqtk seq -A ../../nextflow/data/passenger/test_run/unique_with_at_least_two_reads/2B_S2_001.fastq > 2B_S2_001.fasta
# makeblastdb -in 2B_S2_001.fasta -dbtype nucl -out 2B_S2_001.db
# blastn -max_target_seqs 100000000 -query around-deleted-corrected.fasta -db 2B_S2_001.db -task "blastn-short" -outfmt 6 -evalue 1000 -word_size 7 -perc_identity 95 > 2B_S2_001.blast.tsv

for file in ../../nextflow/data/passenger/test_run/unique_with_at_least_two_reads/*.fastq
do
    filename=$(basename "$file" .fastq)
    seqtk seq -A "$file" > "$filename.fasta"
    makeblastdb -in "$filename.fasta" -dbtype nucl -out "$filename.db"
    blastn -max_target_seqs 100000000 -query around-deleted-corrected.fasta -db "$filename.db" -task "blastn-short" -outfmt 6 -evalue 1000 -word_size 7 -perc_identity 95 > "$filename.blast.tsv"
done
