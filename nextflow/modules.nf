nextflow.enable.dsl =2

/*
 * Process 1A: trim the first three bases of the paired end reads.
 */

// meh, let's un-capitalize these process publishDir's
process TRIM_COMBINE_MATES { 
  time '30m'
  memory '2g'
  cpus 8
  cache 'lenient'
  container 'quay.io/matsengrp/gcreplay-pipeline:latest'
  publishDir "$params.results/TRIM_COMBINE_MATES/" 
  
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
// TODO why are we including the date in the prefix? Seems unnecessary
process DEMULTIPLEX_PLATES {
  time '10m'
  memory '2g'
  cpus 1
  cache 'lenient'
  container 'quay.io/matsengrp/gcreplay-pipeline:latest'
  publishDir "$params.results/DEMULTIPLEX_PLATES/" 

  input: 
    tuple val(key), val(key_file), val(date), path(key_fasta)
    path(plate_barcodes)
  output: tuple val(key), val(key_file), path("${key}.${date}.*")
  script:
  """
  cat ${key_fasta} | fastx_barcode_splitter.pl \
    --bcfile ${plate_barcodes} --eol \
    --prefix ${key}.${date}. --exact
  """
}

/*
 * Process 1C: Demultiplex each plate into the 96 wells per plate
 */
process DEMULTIPLEX_WELLS {
  time '10m'
  memory '2g'
  cpus 1
  cache 'lenient'
  container 'quay.io/matsengrp/gcreplay-pipeline:latest'
  publishDir "$params.results/DEMULTIPLEX_WELLS/"

  input: 
    tuple val(key), val(key_file), path(plate)
    path(well_barcodes)
  output: tuple val(key), val(key_file), path("${plate}.*")
  script:
  """
  cat ${plate} | fastx_barcode_splitter.pl \
    --bcfile ${well_barcodes} --bol --prefix ${plate}.
  """
}
// --bcfile ${params.reads_prefix}/${params.well_barcodes} --bol --prefix ${plate}.


/*
 * Process 1D: Split each demultiplexed fasta into 
 * heavy and light chain by using cutadapt to search
 * for common motifs
 */
process SPLIT_HK {
  time '20m'
  memory '2g'
  cpus 4
  cache 'lenient'
  container 'quay.io/matsengrp/gcreplay-pipeline:latest'
  publishDir "$params.results/SPLIT_HK/"

  input: 
    tuple val(key), val(key_file), path(well) 
    val(motif)
    val(chain)
  output: tuple val(key), val(key_file), path("${well}.${chain}")
  script:
  """
  cutadapt --cores ${task.cpus} -g ${motif} -e 0.2 ${well} --discard-untrimmed -o ${well}.${chain}
  """
}


/*
 * Process 1E: Collapse each heavy and light chain
 * demultiplexed files into  
 */
process COLLAPSE_RANK_PRUNE {
  time '20m'
  memory '2g'
  cpus 1
  cache 'lenient'
  container 'quay.io/matsengrp/gcreplay-pipeline:latest'
  publishDir "$params.results/COLLAPSE_RANK_PRUNE/"

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
  time '20m'
  memory '2g'
  cpus 1
  cache 'lenient'
  container 'quay.io/matsengrp/gcreplay-pipeline:latest'
  publishDir "$params.results/MERGE_BCRS/"

  input: tuple val(key), val(key_file), path(all_coll_rank)
  output: tuple val(key), val(key_file), path("${key}.fasta")
  script:
  """
  awk '/>/{sub(">","&"FILENAME".")}1' ${all_coll_rank} > merged.fasta
  gcreplay-tools-post-collapse.py sort-fasta --fasta merged.fasta -o ${key}.fasta
  """
}


/*
 * Process 2A: Annotate the top ranked seqs
 */
process PARTIS_ANNOTATION {
  time '30m'
  memory '4g'
  cpus 4
  cache 'lenient'
  container 'quay.io/matsengrp/gcreplay-pipeline:partis'
  publishDir "$params.results/PARTIS_ANNOTATION/"

  input: 
    tuple val(key), val(key_file), path(merged_fasta)
  output: tuple val(key), val(key_file), path(merged_fasta), path("${key}/")
  script:
  """
  wd=\$PWD
  cd /partis
  initial-annotate.sh \${wd}/${merged_fasta} \${wd}/${key} $params.partis_anno_dir
  """
}


/*
 * Process 2B: Wrangle and parse the annotations
 */
process PARTIS_WRANGLE {
  time '15m'
  memory '2g'
  cpus 1
  cache 'lenient'
  container 'quay.io/matsengrp/gcreplay-pipeline:latest'
  publishDir "$params.results/PARTIS_WRANGLE/"

  input: tuple val(key), path(key_file), path(merged_fasta), path(partis_out)
  output: path "annotated-${key}*.csv"
  
  """  
  IGH_AIRR=${partis_out}/engrd/single-chain/partition-igh.tsv
  IGK_AIRR=${partis_out}/engrd/single-chain/partition-igk.tsv

  # wrangle annotation -> gc merged dataframe
  gcreplay-tools-post-collapse.py wrangle-annotation \
      --igh-airr \$IGH_AIRR \
      --igk-airr \$IGK_AIRR \
      --input-fasta $merged_fasta \
      --key-file $key_file \
      -o ${key}-gc-df-hk.csv
  
  # now, split the wrangled df into single mouse / gc
  # --sample 10 TODO add this option to pipeline params
  gcreplay-tools-post-collapse.py df-groupby \
      -df ${key}-gc-df-hk.csv \
      -o annotated-${key}
  """
}


/*
 * Process 3A: Wrangle and featurize nodes
 */
process GCTREE {
  time '8h'
  memory '64g'
  cpus 1
  cache 'lenient'
  container 'quay.io/matsengrp/gcreplay-pipeline:latest'
  publishDir "$params.results/GCTREE/", mode: "copy"

  input: 
    path single_mouse_gc_df
    path hdag_sub
    path hdag_mut
    path dms_vscores
    path dms_sites
  output: path("PR*")
  shell:
  template "gctree_infer_featurize.sh"
}


/*
 * Process 3B: Merge all results
 */
// Maybe rename to BCR_AND_NODE_TABLES
process MERGE_RESULTS {
  time '10m'
  memory '2g'
  cpus 1
  cache 'lenient'
  container 'quay.io/matsengrp/gcreplay-pipeline:latest'
  publishDir "$params.results/MERGE_RESULTS/", mode: "copy"

  input: path(all_results)
  output: tuple path("observed-seqs.csv"), path("gctree-node-data.csv")
  script:
  """
  gcreplay-tools-post-collapse.py merge-results
  """
}


/*
 * Process 4A: NDS-LB analysis notebook papermill
 */
process NDS_LB_ANALYSIS {
  time '10m'
  memory '2g'
  cpus 1
  cache 'lenient'
  container 'quay.io/matsengrp/gcreplay-pipeline:56_analysis_in_pipeline'
  publishDir "$params.results/NDS_LB_ANALYSIS/", mode: "copy"

  input: 
    tuple path(all_results), path(metadata), val(ranking_coeff_subdir), val(svg_scale) 

  output: 
    tuple path("NDS_LB.ipynb"), path("*.pdf"), path("*.csv"), path("*.svg")

  script:
  """
  #! /bin/bash
  activate_env replay
  papermill NDS_LB.ipynb NDS_LB.ipynb \
    -p results './' \
    -p ranking_subdir ${ranking_coeff_subdir} \
    -p scale ${svg_scale} \
    -p metadata_csv ${metadata} \
    -p outbase './'
  """
}



 
