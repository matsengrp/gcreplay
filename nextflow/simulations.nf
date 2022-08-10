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

params.simulations      = "$projectDir/data/simulations"
params.results          = "$projectDir/results"
params.hdag_sub         = "data/mutability/MK_RS5NF_substitution.csv"
params.hdag_mut         = "data/mutability/MK_RS5NF_mutability.csv"
params.dms_vscores      = "https://media.githubusercontent.com/media/jbloomlab/Ab-CGGnaive_DMS/main/results/final_variant_scores/final_variant_scores.csv"
params.dms_sites        = "https://raw.githubusercontent.com/jbloomlab/Ab-CGGnaive_DMS/main/data/CGGnaive_sites.csv"
params.dms_mvscores     = "data/torchdms/prepped-dms/final-multi-variant-scores.csv"
params.tdms_model       = "data/torchdms/prepped-dms/tdms-fully-connected.model"
params.tdms_model_lin   = "data/torchdms/prepped-dms/tdms-linear.model"
params.igk_idx          = 336


log.info """\
G C Re - F L O W (simulations)!
Matsen, Victora Labs
Fred Hutchinson CRC, Seattle WA
Rockefeller University, New York NY.
================================
"""


process GCTREE_SIM {
  container 'quay.io/matsengrp/gcreplay-pipeline:latest'
  publishDir "$params.results/simulation-gctrees/", mode: "copy"
  label "mem_large"
  input: path(simulated_gc)
  output: path("gct-")
  shell:
  template "gctree_infer_featurize_simu.sh"
}


workflow {

    simu = Channel.fromPath("$params.simulations")
    simu | GCTREE_SIM

}




