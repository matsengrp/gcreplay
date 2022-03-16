nextflow.enable.dsl =2

/*
 * Process 1A: trim the first three bases of the paired end reads.
 */
process TRIM_COMBINE_MATES { 
  container 'quay.io/matsengrp/gcreplay-pipeline:2022-03-03'
  publishDir "$params.results/trimmed_combined_fasta/" 
  label 'multithread'
  input: tuple val(key), val(key_file), val(date), path(read1), path(read2)
  output: tuple val(key), val(key_file), val(date), path("${key}.fasta")
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
  container 'quay.io/matsengrp/gcreplay-pipeline:2022-03-03'
  publishDir "$params.results/demultiplexed_plates_fasta/" 
  input: tuple val(key), val(key_file), val(date), path(key_fasta)
  output: tuple val(key), val(key_file), path("${key}.${date}.*")
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
  container 'quay.io/matsengrp/gcreplay-pipeline:2022-03-03'
  publishDir "$params.results/demultiplexed_wells_fasta/"
  input: tuple val(key), val(key_file), path(plate)
  output: tuple val(key), val(key_file), path("${plate}.*")
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
  container 'quay.io/matsengrp/gcreplay-pipeline:2022-03-03'
  publishDir "$params.results/split_HK/"
  label 'multithread'
  input: 
    tuple val(key), val(key_file), path(well) 
    val(motif)
    val(chain)
  output: tuple val(key), val(key_file), path("${well}.${chain}")
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
  container 'quay.io/matsengrp/gcreplay-pipeline:2022-03-03'
  publishDir "$params.results/rank_collapsed/"
  input: tuple val(key), val(key_file), path(well_chain)
  output: tuple val(key), val(key_file), path("${well_chain}.R")
  script:
  """
  fastx_collapser -i ${well_chain} -o rank_collapsed.fasta
  gcreplay-tools.py curate-high-count-seqs --fasta rank_collapsed.fasta \
    --count-threshold ${params.bcr_count_thresh} \
    -o ${well_chain}.R
  """
}


/*
 * Process 1F: Merge the top ranked BCR's
 */
process MERGE_BCRS {
  container 'quay.io/matsengrp/gcreplay-pipeline:2022-03-03'
  publishDir "$params.results/ranked_bcr_sequences_per_well/"
  input: tuple val(key), val(key_file), path(all_coll_rank)
  output: tuple val(key), val(key_file), path("${key}.fasta")
  script:
  """
  awk '/>/{sub(">","&"FILENAME".")}1' ${all_coll_rank} > merged.fasta
  gcreplay-tools.py sort-fasta --fasta merged.fasta -o ${key}.fasta
  """
}


/*
 * Process 2A: Annotate the top ranked seqs
 */
process PARTIS_ANNOTATION {
  container 'quay.io/matsengrp/partis:main'
  publishDir "$params.results/partis_annotation/"
  input: tuple val(key), val(key_file), path(merged_fasta)
  output: tuple val(key), val(key_file), path(merged_fasta), path("${key}/")
  script:
  """
  wd=\$PWD
  cd /partis
  initial-annotate.sh \${wd}/${merged_fasta} \${wd}/${key} ${params.partis_anno_dir}germlines/
  """
}


/*
 * Process 2B: Wrangle and parse the annotations
 */
process PARTIS_WRANGLE {
  container 'quay.io/matsengrp/gcreplay-pipeline:2022-03-03'
  publishDir "$params.results/single_gc_wrangle/"
  input: tuple val(key), path(key_file), path(merged_fasta), path(partis_out)
  output: path "annotated-${key}*.csv"
  
  """  
  IGH_AIRR=${partis_out}/engrd/single-chain/partition-igh.tsv
  IGK_AIRR=${partis_out}/engrd/single-chain/partition-igk.tsv

  # wrangle annotation -> gc merged dataframe
  gcreplay-tools.py wrangle-annotation \
      --igh-airr \$IGH_AIRR \
      --igk-airr \$IGK_AIRR \
      --input-fasta $merged_fasta \
      --key-file $key_file \
      -o ${key}-gc-df-hk.csv
  
  # now, split the wrangled df into single mouse / gc
  gcreplay-tools.py df-groupby \
      -df ${key}-gc-df-hk.csv \
      -o annotated-${key}
  """
  //--sample 10 \
}


/*
 * Process 3A: Wrangle and featurize nodes
 */
process GCTREE {
  container 'quay.io/matsengrp/gcreplay-pipeline:2022-03-03'
  publishDir "$params.results/gctrees/", mode: "copy"
  label "mem_large"
  //errorStrategy 'ignore'
  input: path(single_mouse_gc_df)
  output: path("PR*")
  shell:
  template "gctree_infer_featurize.sh"
}


/*
 * Process 3B: Merge all results
 */
process MERGE_RESULTS {
  container 'quay.io/matsengrp/gcreplay-pipeline:2022-03-03'
  publishDir "$params.results/merged-results/", mode: "copy"
  label "mem_large"
  input: path(all_results)
  output: tuple path("observed-seqs.csv"), path("gctree-node-data.csv")
  shell:
  """
  gcreplay-tools.py merge-results
  """
}






