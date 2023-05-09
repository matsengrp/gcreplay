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

params.seqs    = "$projectDir/data/beast/test-observed-seqs.fasta"
params.beast_template   = "$projectDir/data/beast/beast_templates/skyline_histlog.template"
params.dms_vscores      = "https://media.githubusercontent.com/media/jbloomlab/Ab-CGGnaive_DMS/main/results/final_variant_scores/final_variant_scores.csv"
params.dms_sites        = "https://raw.githubusercontent.com/jbloomlab/Ab-CGGnaive_DMS/main/data/CGGnaive_sites.csv"
params.results          = "$projectDir/results"
params.chain_length     = 50000000
params.log_every        = 10000
params.burn_frac        = 0.9

log.info """\
G C Re - F L O W (beast)!
Matsen, Victora Labs
Fred Hutchinson CRC, Seattle WA
Rockefeller University, New York NY.
================================
"""


process ADD_TIME_TO_FASTA {
  container 'quay.io/matsengrp/gcreplay-pipeline:beagle-beast-2023-04-24'
  publishDir "$params.results/beast-input/", mode: "copy"
  input: path(seqs)
  output: path("*_with_time.fasta")
  shell:
  """
  add_date_to_fasta.py $seqs
  """
}

process BEAST_TIMETREE {
  label 'large'
  stageInMode 'copy' // I guess beast doesn't like symlinks
  container 'quay.io/matsengrp/gcreplay-pipeline:beagle-beast-2023-04-24'
  publishDir "$params.results/beast"
  input: tuple path(seqs_with_time), path(beast_template)
  output: path("btt-*")
  shell:
  template "beast_time_tree.sh"
}

process ETE_CONVERSION {
  errorStrategy 'retry'
  maxRetries 3
  label 'mem_large'
  container 'quay.io/matsengrp/gcreplay-pipeline:historydag-ete-2023-04-24'
  publishDir "$params.results/ete", mode: "copy" 
  input: path(beast_output)
  output: path("ete-*")
  shell:
  """
  BEAST_OUTDIR=$beast_output
  OUTDIR=ete-\${BEAST_OUTDIR#"btt-"}
  beast2ete.py \
    --xml_file $beast_output/beastgen.xml \
    --nexus_file $beast_output/*.history.trees \
    --dms_df $params.dms_vscores \
    --pos_df $params.dms_sites \
    --burn_frac $params.burn_frac \
    --outdir \$OUTDIR
  """
}


workflow {

    obs = Channel.fromPath("$params.seqs")
    obs | ADD_TIME_TO_FASTA | set {obs_w_time}
    beast_template = Channel.fromPath("$params.beast_template")
    obs_w_time.combine(beast_template) | BEAST_TIMETREE 
    BEAST_TIMETREE.out | ETE_CONVERSION

}




