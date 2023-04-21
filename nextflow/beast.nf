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

params.observed_seqs    = "$projectDir/data/beast/test-observed-seqs.fasta"
// params.observed_seqs    = "$projectDir/data/beast/*.fasta"
params.beast_template   = "$projectDir/data/beast/beast_templates/exponential_histlog.template"
params.results          = "$projectDir/results"

log.info """\
G C Re - F L O W (beast)!
Matsen, Victora Labs
Fred Hutchinson CRC, Seattle WA
Rockefeller University, New York NY.
================================
"""


process ADD_TIME_TO_FASTA {
  container 'e3286480f1a3'
  publishDir "$params.results/beast-timetrees/", mode: "copy"
  // label "mem_large"
  input: path(observed_seqs)
  output: path("*_with_time.fasta")
  shell:
  """
  add_date_to_fasta.py $observed_seqs
  """
}

process BEAST_TIMETREE {
  stageInMode 'copy'
  container 'e3286480f1a3'
  publishDir "$params.results/beast-timetrees/", mode: "copy"
  input: tuple path(observed_seqs_with_time), path(beast_template)
  output: path("btt-*")
  shell:
  template "beast_time_tree.sh"
}


workflow {

    obs = Channel.fromPath("$params.observed_seqs")
    obs | ADD_TIME_TO_FASTA | set {obs_w_time}
    beast_template = Channel.fromPath("$params.beast_template")
    obs_w_time.combine(beast_template) | BEAST_TIMETREE 

}




