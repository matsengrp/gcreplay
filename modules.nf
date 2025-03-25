nextflow.enable.dsl = 2

/*
 * trim the first three bases of the 
 * paired end reads, and combine the mates
 * using pandaseq.
 */
process TRIM_COMBINE_MATES { 
  time '1h'
  memory '8g'
  cpus 16
  cache 'lenient'
  container 'quay.io/matsengrp/gcreplay-pipeline:latest'
  
  input: tuple val(ngs_id), val(date), path(read1), path(read2)
  output: tuple val(ngs_id), val(date), path("${ngs_id}.fasta")
  script:
  """
  gunzip -c ${read1} > ${ngs_id}_R1.fastq
  gunzip -c ${read2} > ${ngs_id}_R2.fastq
  fastx_trimmer -Q33 -i ${ngs_id}_R1.fastq -f 3 -o t_${ngs_id}_R1.fastq
  fastx_trimmer -Q33 -i ${ngs_id}_R2.fastq -f 3 -o t_${ngs_id}_R2.fastq
  pandaseq -T ${task.cpus} -f t_${ngs_id}_R1.fastq -r t_${ngs_id}_R2.fastq -O 0 -w ${ngs_id}.fasta
  """
}


/*
 * demultiplex by plate barcodes.
 * This process takes the combined fasta file
 * from the previous step, and splits it into
 * separate files for each plate.
 */
process DEMULTIPLEX_PLATES {
  time '10m'
  memory '2g'
  cpus 1
  cache 'lenient'
  container 'quay.io/matsengrp/gcreplay-pipeline:latest'

  input: 
    tuple val(ngs_id), val(date), path(ngs_id_fasta)
    path(plate_barcodes)
  output: tuple val(ngs_id), path("${ngs_id}.${date}.*")
  script:
  """
  cat ${ngs_id_fasta} | fastx_barcode_splitter.pl \
    --bcfile ${plate_barcodes} \
    --eol \
    --prefix \
    ${ngs_id}.${date}. \
    --exact
  """
}

/*
 * Demultiplex each plate into the 96 wells per plate
 */
process DEMULTIPLEX_WELLS {
  time '10m'
  memory '2g'
  cpus 1
  cache 'lenient'
  container 'quay.io/matsengrp/gcreplay-pipeline:latest'
  // publishDir "$params.results/DEMULTIPLEX_WELLS/"

  input: 
    tuple val(ngs_id), path(plate)
    path(well_barcodes)
  output: tuple val(ngs_id), path("${plate}.*")
  script:
  """
  cat ${plate} | fastx_barcode_splitter.pl \
    --bcfile ${well_barcodes} \
    --bol \
    --prefix ${plate}.
  """
}


/*
 * Split each demultiplexed fasta into 
 * heavy and light chain by using cutadapt to search
 * for common motifs.
 */
process SPLIT_HK {
  time '20m'
  memory '2g'
  cpus 4
  cache 'lenient'
  container 'quay.io/matsengrp/gcreplay-pipeline:latest'
  // publishDir "$params.results/SPLIT_HK/"

  input: 
    tuple val(ngs_id), path(well) 
    val(motif)
    val(chain)
  output: tuple val(ngs_id), path("${well}.${chain}")
  script:
  """
  cutadapt \
    --cores ${task.cpus} \
    -g ${motif} \
    -e 0.2 ${well} \
    --discard-untrimmed \
    -o ${well}.${chain}
  """
}


/*
 * Collapse each heavy and light chain
 * demultiplexed file, counting the number of occurrences
 * of each unique sequence. This process will also
 * prune low abundance BCR's based on the user
 * defined threshold (default 5).   
 */
process COLLAPSE_RANK_PRUNE {
  time '20m'
  memory '8g'
  cpus 1
  cache 'lenient'
  container 'quay.io/matsengrp/gcreplay-pipeline:latest'
  publishDir "$params.results/COLLAPSE_RANK_PRUNE/"

  input: tuple val(ngs_id), path(well_chain)
  output: tuple val(ngs_id), path("${well_chain}.R")
  script:
  """
  fastx_collapser -i ${well_chain} -o rank_collapsed.fasta
  prune-low-abundance-bcrs.py \
    --fasta rank_collapsed.fasta \
    --count-threshold ${params.bcr_count_thresh} \
    -o ${well_chain}.R
  """
}


/*
 * Merge the top ranked BCR's
 * from all wells into a single fasta file.
 */
process MERGE_BCRS {
  time '20m'
  memory '2g'
  cpus 1
  cache 'lenient'
  container 'quay.io/matsengrp/gcreplay-pipeline:latest'
  publishDir "$params.results/MERGE_BCRS/"

  input: tuple val(ngs_id), path(all_coll_rank)
  output: tuple val(ngs_id), path("${ngs_id}.fasta")
  script:
  """
  awk '/>/{sub(">","&"FILENAME".")}1' ${all_coll_rank} > merged.fasta
  sort-fasta.py --fasta merged.fasta -o ${ngs_id}.fasta
  """
}


/*
 * Annotate the top ranked seqs using partis.
 * This process will take the merged fasta file
 */
process PARTIS_ANNOTATION {
  time '30m'
  memory '4g'
  cpus 4
  cache 'lenient'
  container 'quay.io/matsengrp/gcreplay-pipeline:partis'
  // publishDir "$params.results/PARTIS_ANNOTATION/"

  input: 
    tuple val(ngs_id), path(merged_fasta)
    path(partis_anno_dir)
  output: tuple val(ngs_id), path(merged_fasta), path("${ngs_id}/")
  script:
  """
  initial-annotate.sh $merged_fasta $ngs_id $partis_anno_dir
  """
}


/*
 * Wrangle and parse the annotations,
 * merging the correcsponding well-specific 
 * heavy and light chain
 * sequences into a single dataframe per GC.
 */
process PARTIS_WRANGLE {
  time '15m'
  memory '2g'
  cpus 1
  cache 'lenient'
  container 'quay.io/matsengrp/gcreplay-pipeline:latest'
  // publishDir "$params.results/PARTIS_WRANGLE/"

  input: 
    tuple val(ngs_id), path(merged_fasta), path(partis_out)
    path(gc_metadata)
  output: path "*.csv"
  
  """  
  IGH_AIRR=${partis_out}/engrd/single-chain/partition-igh.tsv
  IGK_AIRR=${partis_out}/engrd/single-chain/partition-igk.tsv

  # wrangle annotation -> gc merged dataframe
  wrangle-annotation.py \
      --igh-airr \$IGH_AIRR \
      --igk-airr \$IGK_AIRR \
      --input-fasta $merged_fasta \
      --gc-metadata $gc_metadata \
      --ngs-id $ngs_id \
      -o gc-df-hk.csv
  
  # now, split the wrangled df into single mouse / gc
  df-groupby.py \
      --dataframe gc-df-hk.csv \
      --column "uid"

  rm gc-df-hk.csv
  """
}

/*
 * Wrangle and featurize nodes
 */
process GCTREE {
  time '8h'
  memory '64g'
  cpus 1
  cache 'lenient'
  container 'quay.io/matsengrp/gcreplay-pipeline:latest'
  publishDir "$params.results/gctrees/", mode: "copy"

  input: 
    path single_mouse_gc_df
    path hdag_sub
    path hdag_mut
    path dms_vscores
    path dms_sites
    path gctree_tools
    path trees_utils
  output: path("${single_mouse_gc_df.baseName}"), type: 'dir'
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
  publishDir "$params.results/", mode: "copy"

  input: path(all_results)
  output: tuple path("observed-seqs.csv"), path("gctree-node-data.csv")
  script:
  """
  merge-results.py
  """
}


/*
 * Process 4A: NDS-LB analysis notebook papermill
 */
process NDS_LB_ANALYSIS {
  time '10m'
  memory '16g'
  cpus 4
  // cache 'lenient'
  container 'quay.io/matsengrp/gcreplay-pipeline:analysis-notebooks'
  publishDir "$params.results/notebooks/NDS-LB/", mode: "copy"

  input: 
    path notebook
    path utils
    path metadata
    tuple path(gctree_dirs), val(ranking_coeff_subdir) //, val(svg_scale) 

  output: 
    tuple(
      val(ranking_coeff_subdir),
      path("$ranking_coeff_subdir/data.csv"),
      path("$ranking_coeff_subdir/$notebook"), 
      path("$ranking_coeff_subdir/scatter*.svg"), 
      path("$ranking_coeff_subdir/stacked_trees/*.svg"), 
      path("$ranking_coeff_subdir/*.html")
    )

  script:
  """
  # move all gctree directories to a single directory
  # to organize, and allow notebooks to be still use the ut.trees module
  mkdir -p gctrees/
  for dir in ${gctree_dirs}; do
    mv \$dir gctrees/
  done

  mkdir -p $ranking_coeff_subdir

  # run the notebook
  papermill $notebook $ranking_coeff_subdir/$notebook \
    -p results '.' \
    -p ranking_subdir $ranking_coeff_subdir \
    -p metadata_csv $metadata \
    -p outbase '.' \
    -p workflow_env_exec True
  """
}
//-p scale $svg_scale \


/*
 * Process 4B: fitness-regression analysis notebook papermill
 */
process FITNESS_REGRESSION_ANALYSIS {
  time '10m'
  memory '16g'
  cpus 1
  // cache 'lenient'
  container 'quay.io/matsengrp/gcreplay-pipeline:analysis-notebooks'
  publishDir "$params.results/notebooks/fitness-regression/", mode: "copy"

  input: 
    path notebook
    path utils
    path metadata
    tuple path(gctree_dirs), val(ranking_coeff_subdir)

  output: 
    tuple(
      path("$ranking_coeff_subdir/$notebook"), 
      path("$ranking_coeff_subdir/*.pdf"),
      path("$ranking_coeff_subdir/*.csv")
    )

  script:
  """
  mkdir -p gctrees/
  for dir in ${gctree_dirs}; do
    mv \$dir gctrees/
  done

  mkdir -p $ranking_coeff_subdir

  # run the notebook
  papermill $notebook $ranking_coeff_subdir/$notebook \
    -p results '.' \
    -p ranking_subdir $ranking_coeff_subdir \
    -p metadata_csv $metadata \
    -p outbase '.'
  """
}

/*
 * Process 4C: mutations analysis notebook papermill
 */
process MUTATIONS_ANALYSIS {
  time '30m'
  memory '16g'
  cpus 4
  // cache 'lenient'
  container 'quay.io/matsengrp/gcreplay-pipeline:analysis-notebooks'
  publishDir "$params.results/notebooks/mutations/", mode: "copy"

  input: 
    path notebook
    path utils
    path metadata
    path dms_vscores
    path dms_sites
    path chigy_hc_mut_rates
    path chigy_lc_mut_rates
    path pdb
    tuple path(gctree_dirs), val(ranking_coeff_subdir)

  output: 
    tuple(
      path("$ranking_coeff_subdir/$notebook"), 
      path("$ranking_coeff_subdir/*.pdf"),
      path("$ranking_coeff_subdir/data.csv"),
      val(ranking_coeff_subdir)
    )

  script:
  """
  mkdir -p gctrees/
  for dir in ${gctree_dirs}; do
    mv \$dir gctrees/
  done

  mkdir -p $ranking_coeff_subdir

  # run the notebook
  papermill $notebook $ranking_coeff_subdir/$notebook \
    -p final_variant_scores $dms_vscores \
    -p dms_sites $dms_sites \
    -p chigy_hc_mut_rates $chigy_hc_mut_rates \
    -p chigy_lc_mut_rates $chigy_lc_mut_rates \
    -p pdb $pdb \
    -p results '.' \
    -p ranking_subdir $ranking_coeff_subdir \
    -p metadata_csv $metadata \
    -p outbase '.' \
    -p workflow_env_exec True
  """
}

/*
 * Process 4D: interactive heatmaps notebook
 */
process INTERACTIVE_HEATMAPS {
  time '10m'
  memory '16g'
  cpus 4
  // cache 'lenient'
  container 'quay.io/matsengrp/gcreplay-pipeline:heatmaps-notebook'
  publishDir "$params.results/notebooks/interactive-heatmaps/", mode: "copy"

  input: 
    path notebook
    tuple path(mutations_data), val(ranking_coeff_subdir)

  output: 
    tuple(
      path("$ranking_coeff_subdir/$notebook"), 
      path("$ranking_coeff_subdir/*.html")
    )

  script:
  """
  mkdir -p $ranking_coeff_subdir

  # run the notebook
  papermill $notebook $ranking_coeff_subdir/$notebook \
    -p mutations_data $mutations_data \
    -p outbase '$ranking_coeff_subdir'
  """
}


/*
 * Process 4E: cell summaries notebook
 */
process CELL_SUMMARIES {
  time '10m'
  memory '16g'
  cpus 1
  // cache 'lenient'
  // stageInMode 'copy'
  container 'quay.io/matsengrp/gcreplay-pipeline:analysis-notebooks'
  publishDir "$params.results/notebooks/cell-summaries/", mode: "copy"

  input: 
    path notebook
    path metadata
    tuple val(ranking_coeff_subdir), path(nds_lb_summary_csv), path(observed_seqs)

  output: 
    tuple(
      path("$ranking_coeff_subdir/$notebook"), 
      path("$ranking_coeff_subdir/*.pdf"),
      path("$ranking_coeff_subdir/*.csv")
    )

  script:
  """
  mkdir -p $ranking_coeff_subdir
  
  # export IPYTHONDIR=/tmp/.ipython
  # mkdir -p \$IPYTHONDIR

  papermill $notebook $ranking_coeff_subdir/$notebook \
    -p metadata_csv $metadata \
    -p cell_table_csv $observed_seqs \
    -p ranking_subdir $ranking_coeff_subdir \
    -p nds_lb_summary_csv $nds_lb_summary_csv \
    -p outbase '.' \
    -p workflow_env_exec True
  """
}

/*
 * 
 */
process PHENOTYPE_TRAJECTORIES {
  time '30m'
  memory '16g'
  cpus 4
  container 'quay.io/matsengrp/gcreplay-pipeline:analysis-notebooks'
  publishDir "$params.results/notebooks/phenotype-trajectories/", mode: "copy"

  input: 
    path notebook
    path utils
    path metadata
    path dms_vscores
    path dms_sites
    path hdag_mut
    path hdag_sub
    tuple path(gctree_dirs), val(ranking_coeff_subdir)

  output: 
    tuple(
      path("$ranking_coeff_subdir/$notebook"), 
      path("$ranking_coeff_subdir/*.pdf"),
      path("$ranking_coeff_subdir/*.csv"),
      path("$ranking_coeff_subdir/*.png"),
      path("$ranking_coeff_subdir/phenotype-trajectories/**"),
      val(ranking_coeff_subdir)
    )

  script:
  """
  mkdir -p gctrees/
  for dir in ${gctree_dirs}; do
    mv \$dir gctrees/
  done

  mkdir -p $ranking_coeff_subdir

  # run the notebook
  papermill $notebook $ranking_coeff_subdir/$notebook \
    -p results '.' \
    -p ranking_subdir $ranking_coeff_subdir \
    -p metadata_csv $metadata \
    -p final_variant_scores $dms_vscores \
    -p dms_sites $dms_sites \
    -p mutability_csv $hdag_mut \
    -p substitution_csv $hdag_sub \
    -p outbase '.' \
    -p workflow_env_exec True
  """
}

/*
 * 
 */
process ANALYSIS_10X {
  time '30m'
  memory '16g'
  cpus 4
  stageInMode 'copy'
  container 'quay.io/matsengrp/gcreplay-pipeline:analysis-notebooks'
  publishDir "$params.results/notebooks/10x/", mode: "copy"

  input: 
    path notebook
    path data_AV1, stageAs: 'data_AV1'
    path data_AV2, stageAs: 'data_AV2'
    path data_AV3, stageAs: 'data_AV3'
    path data_10wk, stageAs: 'data_10wk'
    path dms_vscores
    path dms_sites
    path metadata_sheet
    path metadata_10wk_sheet

  output: 
    tuple(
      path("$notebook"), 
      path("*.pdf"),
      path("*.csv")
    )

  script:
  """
  # run the notebook
  papermill $notebook $notebook \
    -p data_AV1 $data_AV1 \
    -p data_AV2 $data_AV2 \
    -p data_AV3 $data_AV3 \
    -p data_10wk $data_10wk \
    -p final_variant_scores $dms_vscores \
    -p dms_sites $dms_sites \
    -p metadata_sheet $metadata_sheet \
    -p metadata_10wk_sheet $metadata_10wk_sheet \
    -p outbase '.' \
  """
}