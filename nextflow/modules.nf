nextflow.enable.dsl =2

/*
 * Process 1A: trim the first three bases of the paired end reads.
 */
process TRIM_COMBINE_MATES { 
  container 'quay.io/matsengrp/gcreplay-pipeline:trim_combine_demultiplex'
  publishDir "$params.results/trimmed_combined_fasta/" 
  label 'multithread'
  input: tuple val(key), val(date), path(read1), path(read2)
  output: tuple val(key), val(date), path("${key}.fasta")
  script:
  """
  fastx_trimmer -Q33 -i ${read1} -f 3 -o t_${read1}
  fastx_trimmer -Q33 -i ${read2} -f 3 -o t_${read2}
  pandaseq -T 8 -f t_${read1} -r t_${read2} -O 0 -w ${key}.fasta
  """
}


/*
 * Process 1B: demultiplex the different plates
 * if there are no sequences matching a barcode,
 * and thus th resulting file is empty, we don't
 * output it.
 */
//path plate_barcodes
process DEMULTIPLEX_PLATES {
  container 'quay.io/matsengrp/gcreplay-pipeline:trim_combine_demultiplex'
  publishDir "$params.results/demultiplexed_plates_fasta/" 
  input: tuple val(key), val(date), path(key_fasta)
  output: tuple val(key), path("${key}.${date}.*")
  script:
  """
  cat ${key_fasta} | fastx_barcode_splitter.pl \
    --bcfile ${params.reads_prefix}/${params.plate_barcodes} --eol \
    --prefix ${key}.${date}. --exact
  """
}

/*
 * Process 1C: Demultiplex each plate into the 96 wells per plate
 */
process DEMULTIPLEX_WELLS {
  container 'quay.io/matsengrp/gcreplay-pipeline:trim_combine_demultiplex'
  publishDir "$params.results/demultiplexed_wells_fasta/"
  input: tuple val(key), path(plate)
  output: tuple val(key), path("${plate}.*")
  script:
  """
  cat ${plate} | fastx_barcode_splitter.pl \
    --bcfile ${params.reads_prefix}/${params.well_barcodes} --bol --prefix ${plate}.
  """
}


/*
 * Process 1D: Split each demultiplexed fasta into 
 * heavy and light chain by using cutadapt to search
 * for common motifs
 */
process SPLIT_HK {
  container 'quay.io/matsengrp/gcreplay-pipeline:trim_combine_demultiplex'
  publishDir "$params.results/split_HK/"
  label 'multithread'
  input: 
    tuple val(key), path(well) 
    val(motif)
    val(chain)
  output: tuple val(key), path("${well}.${chain}")
  script:
  """
  echo \$PATH
  cutadapt --cores 8 -g ${motif} -e 0.2 ${well} --discard-untrimmed -o ${well}.${chain}
  """
}


/*
 * Process 1E: Collapse each heavy and light chain
 * demultiplexed files into  
 */
process COLLAPSE_RANK_PRUNE {
  container 'quay.io/matsengrp/gcreplay-pipeline:trim_combine_demultiplex'
  publishDir "$params.results/rank_collapsed/"
  input: tuple val(key), path(well_chain)
  output: tuple val(key), path("${well_chain}.R")
  script:
  if( params.top_n_rank != 0 )
    """
    fastx_collapser -i ${well_chain} | head -n ${params.top_n_rank} > ${well_chain}.R
    """
}

// TODO ^^  you dont really need to do this (below) .. until we care about getting thrown
// Sequences
// 
//else if( params.n_threshold != 0 )
//  """
//  fastx_collapser -i ${well_chain} > collapsed_well_chain.fasta
//  ./threshold-counts.py --fasta collapsed_well_chain.fasta --out ${well_chain}.R
//  """


/*
 * Process 1F: Merge the top ranked BCR's
 */
process MERGE_BCRS {
  container 'quay.io/matsengrp/gcreplay-pipeline:trim_combine_demultiplex'
  publishDir "$params.results/ranked_bcr_sequences_per_well/"
  input: tuple val(key), path(all_coll_rank)
  output: tuple val(key), path("${key}.fasta")
  script:
  """
  awk '/>/{sub(">","&"FILENAME".")}1' ${all_coll_rank} > ${key}.fasta
  """
}


/*
 * Process 2A: Annotate the top ranked seqs
 */
process PARTIS_ANNOTATION {
  container 'quay.io/matsengrp/partis:dev'
  publishDir "$params.results/partis_annotation/"
  input: tuple val(key), path(merged_fasta)
  output: path("${key}/")
  script:
  """
  wd=\$PWD
  cd /partis
  initial-annotate.sh \${wd}/${merged_fasta} \${wd}/${key} ${params.partis_anno_dir}germlines/
  """
}
//cd /partis && mkdir \${wd}/${key}
