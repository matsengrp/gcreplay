/*
 * This Source Code Form is subject to the terms of the GNU GENERAL PUBLIC LICENCE
 * License, v. 3.0.
 */


/*
 * 'beast-flow' Run beast on gcreplay data
 * 
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
params.beast_template   = "$projectDir/data/beast/beast_templates/skyline_histlog.template.patch"
params.results          = "$projectDir/results"
params.chain_length     = 50000000
params.log_every        = 10000
params.convert_to_ete   = true
params.dms_vscores      = "https://media.githubusercontent.com/media/jbloomlab/Ab-CGGnaive_DMS/improved-Kd-fitting/tite-seq-modeling/output/final_variant_scores.csv"
params.dms_sites        = "https://raw.githubusercontent.com/jbloomlab/Ab-CGGnaive_DMS/main/data/CGGnaive_sites.csv"
params.burn_frac        = 0.9
params.save_pkl_trees   = false

log.info """\
G C Re - F L O W (beast)!
Matsen, Victora Labs
Fred Hutchinson CRC, Seattle WA
Rockefeller University, New York NY.
================================
"""

// TODO add naive sequence to the fasta's?

process ADD_TIME_TO_FASTA {
  time '5m'
  memory '2g'
  cpus 1
  container 'quay.io/matsengrp/gcreplay-pipeline:beagle-beast-2023-04-24'
  publishDir "$params.results/beast-input/", mode: "copy"
  input: path(seqs)
  output: path("*_with_time.fasta")
  shell:
  """
  add_date_to_fasta.py $seqs
  """
}

process BEASTGEN {
  // label 'small'
  time '5m'
  memory '2g'
  cpus 1
  stageInMode 'copy' // I guess beast doesn't like symlinks
  container 'quay.io/matsengrp/gcreplay-pipeline:beagle-beast-2023-04-24'
  input: tuple path(seqs_with_time), path(beast_template)
  output: path("btt-*")
  shell:
  template "beastgen.sh"
}

process NAIVE_ROOT_PATCH {
  // label 'small'
  time '5m'
  memory '2g'
  cpus 1
  container 'quay.io/matsengrp/gcreplay-pipeline:historydag-ete-2023-04-24'
  publishDir "$params.results/beastgen"
  input: path(beastgen_output)
  output: path(beastgen_output)
  shell:
  """
  beast_template_root_fix.py \
    --templated_beast_xml_path $beastgen_output/beastgen.xml \
    --out_path $beastgen_output/beastgen.naiveroot.xml
  """  
}

process BEAST_TIMETREE {
  // label 'threaded'
  time '4h'
  memory '8g'
  cpus 16
  stageInMode 'copy' // I guess beast doesn't like symlinks
  container 'quay.io/matsengrp/gcreplay-pipeline:beagle-beast-2023-04-24'
  publishDir "$params.results/beast", mode: "copy"
  input: path(beastgen_output)
  output: path(beastgen_output)
  shell:
  """
  beast -threads $task.cpus -prefix beast -working $beastgen_output/beastgen.naiveroot.xml
  """

}

process ETE_CONVERSION {
  errorStrategy 'retry'
  maxRetries 1
  label 'mem_large'
  container 'quay.io/matsengrp/gcreplay-pipeline:historydag-ete-2023-04-24'
  publishDir "$params.results/ete"
  input: path(beast_output)
  output: path("ete-*")
  shell:
  """
  BEAST_OUTDIR=$beast_output
  OUTDIR=ete-\${BEAST_OUTDIR#"btt-"}
  beast2ete.py \
    --xml_file $beast_output/beastgen.naiveroot.xml \
    --nexus_file $beast_output/*.history.trees \
    --dms_df $params.dms_vscores \
    --pos_df $params.dms_sites \
    --burn_frac $params.burn_frac \
    --outdir \$OUTDIR \
    --save_pkl_trees $params.save_pkl_trees
  """
}

process MERGE_SLICE_DFS {
  container 'quay.io/matsengrp/gcreplay-pipeline:historydag-ete-2023-04-24'
  input: path(all_ete_outputs)
  publishDir "$params.results/phenotype_trajectory", mode: "copy" 
  output: path("slice_df.csv")
  script:
  """
  #!/usr/bin/env python  
  import glob
  import pandas as pd

  pd.concat(
    [
      pd.read_csv(f"{ete_dir}/slice_df.csv")
      for ete_dir in glob.glob('ete-*')
    ]
  ).to_csv("slice_df.csv", index=False)
  """
}

workflow {

    seqs = Channel.fromPath("$params.seqs", checkIfExists: true)
    beast_template = Channel.fromPath("$params.beast_template", checkIfExists: true)

    seqs | ADD_TIME_TO_FASTA | combine(beast_template) \
        | BEASTGEN 
        | NAIVE_ROOT_PATCH \
        | BEAST_TIMETREE \
        | set {beast_outputs}

    if (params.convert_to_ete)
        ETE_CONVERSION(beast_outputs) \
        | collect | MERGE_SLICE_DFS

}




