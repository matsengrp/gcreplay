/*
 * This Source Code Form is subject to the terms of the GNU GENERAL PUBLIC LICENCE
 * License, v. 3.0.
 */


/*
 * 'gcre-flow' - A Nextflow pipeline for running gc analysis workflow
 * smart-seq analysis
 *
 * Fred Hutchinson Cancer Research Center, Seattle WA.
 * Rockefeller University, New York NY.
 *
 * Jared Galloway
 * Tatsuya Araki
 * Duncan Ralph
 * Will DeWitt
 * Will Dumm
 */



/*
 * Enable DSL 2 syntax
 */
nextflow.enable.dsl = 2


/*
 * Define the default parameters - example data get's run by default
 */

params.reads_prefix     = "$projectDir"
params.hk_df =      "$projectDir/data/smartseq/smartseq_patch.csv"
params.key =        "PR-1-9"
params.hdag_sub         = "data/mutability/MK_RS5NF_substitution.csv"
params.hdag_mut         = "data/mutability/MK_RS5NF_mutability.csv"
params.dms_vscores      = "https://media.githubusercontent.com/media/jbloomlab/Ab-CGGnaive_DMS/main/results/final_variant_scores/final_variant_scores.csv"
params.dms_sites        = "https://raw.githubusercontent.com/jbloomlab/Ab-CGGnaive_DMS/main/data/CGGnaive_sites.csv"
params.igk_idx          = 336


log.info """\
G C Re - F L O W (smart-Seq)!
Matsen, Victora Labs
Fred Hutchinson CRC, Seattle WA
Rockefeller University, New York NY.
================================
"""

/*
 * Process 2B: Wrangle and parse the annotations
 */
process SS_SPLIT_GCS {
  container 'quay.io/matsengrp/gcreplay-pipeline:2022-03-03'
  publishDir "$params.results/single_gc_wrangle/"
  input: 
    val key
    path hk_df
  output: path "annotated-${key}*.csv"
  
  """  
  gcreplay-tools.py df-groupby \
      -df ${hk_df} \
      -o annotated-${key}
  """
}


/*
 * Process 3A: Wrangle and featurize nodes
 */
process SS_GCTREE {
  container 'quay.io/matsengrp/gcreplay-pipeline:2022-03-03'
  publishDir "$params.results/gctrees/", mode: "copy"
  label "mem_large"
  input: path(single_mouse_gc_df)
  output: path("PR*")
  shell:
  template "gctree_infer_featurize.sh"
}


/*
 * Process 3B: Merge all results
 */
process SS_MERGE_RESULTS {
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


workflow {

    hk_df_ch = Channel.fromPath(params.hk_df)
    SS_SPLIT_GCS(params.key, hk_df_ch) \
      | flatten() | SS_GCTREE \
      | collect | SS_MERGE_RESULTS

}




