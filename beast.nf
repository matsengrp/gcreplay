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

params.seqs               = "$projectDir/data/beast/test-observed-seqs.fasta"
params.results            = "$projectDir/results"
params.beast_template     = "$projectDir/data/beast/beast_templates/skyline_histlog.template.patch"
params.dms_vscores        = "$projectDir/data/dms/final_variant_scores.csv"
params.dms_sites          = "$projectDir/data/dms/CGGnaive_sites.csv"
params.chain_length       = 25000000
params.log_every          = 10000
params.convert_to_ete     = true
params.burn_frac          = 0.96
params.save_pkl_trees     = false
params.save_single_pkl    = true
params.naive_seq_time     = 0
params.observed_seq_time  = 1

log.info """\
G C Re - F L O W (beast)!
================================
"""


process ADD_TIME_TO_FASTA {
  time '15m'
  memory '4g'
  cpus 1

  cache 'lenient'  
  container 'quay.io/matsengrp/gcreplay-pipeline:beagle-beast-2023-04-24'
  // publishDir "$params.results/ADD_TIME_TO_FASTA/"
  input: tuple val(id), path(fasta, stageAs: 'input.fasta')
  output: tuple val(id), path("${id}.fasta")
  shell:
  """
  add_date_to_fasta.py \
    --delim "@" \
    --naive_keyword "naive" \
    --naive_seq_time $params.naive_seq_time \
    --observed_seq_time $params.observed_seq_time \
    --output ${id}.fasta \
    input.fasta \
  """
}

process BEASTGEN {
  time '5m'
  memory '2g'
  cpus 1
  
  cache 'lenient'
  stageInMode 'copy' // I guess beast doesn't like symlinks
  container 'quay.io/matsengrp/gcreplay-pipeline:beagle-beast-2023-04-24'
  // publishDir "$params.results/BEASTGEN"
  input: 
    tuple val(id), path(fasta)
    path(beast_template)
  output: tuple val(id), path("${id}.beastgen.xml")

  script:
  """
  beastgen \
    -date_prefix '@' \
    -date_order -1 \
    -D "chain_length=${params.chain_length},log_every=${params.log_every}" \
    $beast_template \
    $fasta \
    ${id}.beastgen.xml
  """

}

process NAIVE_ROOT_PATCH {
  time '5m'
  memory '2g'
  cpus 1

  cache 'lenient' 
  container 'quay.io/matsengrp/gcreplay-pipeline:historydag-ete-2023-04-24'
  // publishDir "$params.results/NAIVE_ROOT_PATCH"
  input: tuple val(id), path(beastgen_xml)
  output: tuple val(id), path("${id}.beastgen.naiveroot.xml")
  shell:
  """
  beast_template_root_fix.py \
    --templated_beast_xml_path $beastgen_xml \
    --out_path ${id}.beastgen.naiveroot.xml
  """  
}

process BEAST_TIMETREE {
  time '4h'
  memory '8g'
  cpus 16

  cache 'lenient' 
  // stageInMode 'copy' // I guess beast doesn't like symlinks
  container 'quay.io/matsengrp/gcreplay-pipeline:beagle-beast-2023-04-24'
  publishDir "$params.results/BEAST_TIMETREE"

  input: tuple val(id), path(beastgen_root_patched)
  output: 
    tuple(
      val(id), 
      path(beastgen_root_patched), 
      path("${id}.history.trees"), 
      path("${id}.trees"),
      path("${id}.log")
    )
  shell:
  """
  beast -threads $task.cpus $beastgen_root_patched
  """

}

process ETE_CONVERSION {  
  time '1h'
  memory '16g'
  cpus 1

  cache 'lenient' 
  container 'quay.io/matsengrp/gcreplay-pipeline:historydag-ete-2023-04-24'
  publishDir "$params.results/ETE_CONVERSION", mode: "copy"

  input: 
    tuple(
      val(id), 
      path(beastgen_root_patched), 
      path(beast_history_trees), 
      path(beast_trees), 
      path(beast_log)
    )
    path(dms_vscores)
    path(dms_sites)

  output: path("${id}/"), type: 'dir'
  shell:
  """
  mkdir -p ${id}
  (
  beast2ete.py \
    --xml_file $beastgen_root_patched \
    --nexus_file $beast_history_trees \
    --dms_df $dms_vscores \
    --pos_df $dms_sites \
    --burn_frac $params.burn_frac \
    --outdir ${id} \
    --save_pkl_trees $params.save_pkl_trees \
    --save_single_pkl $params.save_single_pkl
  ) > ${id}/ete_conversion.log 2>&1
  """
}

process MERGE_SLICE_DFS {
  time '2h'
  memory '16g'
  cpus 1

  cache 'lenient'
  container 'quay.io/matsengrp/gcreplay-pipeline:historydag-ete-2023-04-24'
  publishDir "$params.results/MERGE_SLICE_DFS", mode: "copy" 


  input: path(all_ete_outputs)
  output: path("slice_df.csv.gz")
  script:
  """
  #!/usr/bin/env python  
  import glob
  import pandas as pd

  pd.concat(
    [
      pd.read_csv(f"{ete_dir}/slice_df.csv").assign(
        uid=ete_dir.split("/")[0]
      )
      for ete_dir in glob.glob("*/")
    ]
  )[
    ['uid','tree_idx','entry_idx','time','delta_expr','delta_bind_CGG']
  ].to_csv("slice_df.csv.gz", index=False, compression="gzip")
  """
}

workflow BEAST_FLOW {
    
    take: id_seq_ch

    main:      
      ADD_TIME_TO_FASTA(id_seq_ch) | set{id_seq_w_time_ch}

      BEASTGEN(id_seq_w_time_ch, file("$params.beast_template")) \
        | NAIVE_ROOT_PATCH \
        | BEAST_TIMETREE \
        | set {beast_outputs}
      
      if (params.convert_to_ete)
          ETE_CONVERSION(
            beast_outputs, 
            file("$params.dms_vscores"), 
            file("$params.dms_sites")
          ) \
          | collect | MERGE_SLICE_DFS
    
    // emit:
      // MERGE_SLICE_DFS.out
      // ETE_CONVERSION.out

}

workflow {

    id_seq_ch = Channel.fromPath("$params.seqs", checkIfExists: true) \
      | map{it -> [it.baseName, it]} 

    id_seq_ch | BEAST_FLOW

}
