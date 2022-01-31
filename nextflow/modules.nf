nextflow.enable.dsl =2

/*
 * Process 1A: trim the first three bases of the paired end reads.
 */
process TRIM_COMBINE_MATES { 
  container 'quay.io/matsengrp/gcreplay-pipeline:trim_combine_demultiplex'
  //publishDir 'intermediate/trimmed/'
  label 'multithread'
  input: tuple val(key), path(read1), path(read2)
  output: tuple val(key), path("${key}.fasta")
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
  //publishDir 'intermediate/plate_splits/'
  input: tuple val(key), path(key_fasta)
  output: path "${key}.*"
  script:
  """
  cat ${key_fasta} | fastx_barcode_splitter.pl \
    --bcfile ${params.plate_barcodes} --eol \
    --prefix ${key}. --exact
  """
}

/*
 * Process 1C: Demultiplex each plate into the 96 wells per plate
 */
process DEMULTIPLEX_WELLS {
  container 'quay.io/matsengrp/gcreplay-pipeline:trim_combine_demultiplex'
  //publishDir 'intermediate/well_splits/'
  input: path(plate)
  output: path "${plate}.*"
  script:
  """
  cat ${plate} | fastx_barcode_splitter.pl \
    --bcfile ${params.well_barcodes} --bol --prefix ${plate}. --exact
  """
}


/*
 * Process 1D: Split each demultiplexed fasta into 
 * heavy and light chain by using cutadapt to search
 * for common motifs
 */
process SPLIT_HK {
  container 'quay.io/matsengrp/gcreplay-pipeline:trim_combine_demultiplex'
  //publishDir 'intermediate/hk_splits/'
  label 'multithread'
  input: 
    path(well) 
    val(motif)
    val(chain)
  output: path "${well}.${chain}"
  script:
  """
  cutadapt --cores 8 -g ${motif} -e 0.2 ${well} --discard-untrimmed -o ${well}.${chain}
  """
}


/*
 * Process 1E: Collapse each heavy and light chain
 * demultiplexed files into  
 */
process COLLAPSE_RANK_PRUNE {
  container 'quay.io/matsengrp/gcreplay-pipeline:trim_combine_demultiplex'
  //publishDir 'intermediate/coll_rank/'
  input: path(well_chain)
  output: path("${well_chain}.coll_rank")
  script:
  """
  fastx_collapser -i ${well_chain} | head -n 6 > ${well_chain}.coll_rank
  """
}


/*
 * Process 1F: Merge the top ranked BCR's
 */
process MERGE_BCRS {
  container 'quay.io/matsengrp/gcreplay-pipeline:trim_combine_demultiplex'
  publishDir "$params.results/ranked_bcr_sequences_per_well/"
  input: 
    //val(key)
    path(all_coll_rank)
  //output: path("${key}.fasta")
  output: path("merged.fasta")
  script:
  """
  awk '/>/{sub(">","&"FILENAME"_")}1' ${all_coll_rank} > merged.fasta
  """
}


/*
 * Process 2A: Annotate the top ranked seqs
 */
process PARTIS_ANNOTATION {
  container 'quay.io/matsengrp/partis:dev'
  publishDir 'intermediate/partis_annotation/'
  input: path(merged_fasta)
  output: path("out/engrd/single-chain/*.tsv")
  script:
  """
  wd=\$PWD
  cd /partis
  initial-annotate.sh \${wd}/${merged_fasta} \${wd}/out ${params.partis_anno_dir}germlines/
  """
}
  //cd /partis
// ls bin/partis
// ls /fh/fast/matsen_e/shared/gcreplay/nextflow/data/partis_annotation/germlines/





