/*
 * This Source Code Form is subject to the terms of the GNU GENERAL PUBLIC LICENSE
 * License, v. 3.0.
 */


/*
 * 'gcre-flow' - A Nextflow pipeline for running gc analysis workflow
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

// TODO move to config file
params.ngs_manifest         = "$projectDir/data/input/ngs_manifest.csv"
params.gc_metadata          = "$projectDir/data/input/gc_metadata.csv"
params.reads_prefix         = "$projectDir/data/input/NGS-gz"
params.results              = "$projectDir/results/"

params.plate_barcodes       = "$projectDir/data/barcodes/plateBC.txt"
params.well_barcodes        = "$projectDir/data/barcodes/96FBC.txt"
params.partis_anno_dir      = "$projectDir/data/partis_annotation/germlines"
params.hdag_sub             = "$projectDir/data/mutability/MK_RS5NF_substitution.csv"
params.hdag_mut             = "$projectDir/data/mutability/MK_RS5NF_mutability.csv"
params.chigy_hc_mut_rates   = "$projectDir/data/mutability/chigy_hc_mutation_rates_nt.csv"
params.chigy_lc_mut_rates   = "$projectDir/data/mutability/chigy_lc_mutation_rates_nt.csv"
params.pdb                  = "$projectDir/data/AbCGG_structure/combined_ch2_eh2-coot_IMGT.pdb"
params.dms_vscores          = "$projectDir/data/dms/final_variant_scores.csv"
params.dms_sites            = "$projectDir/data/dms/CGGnaive_sites.csv"
params.heavy_chain_motif    = "aGCgACgGGaGTtCAcagACTGCAACCGGTGTACATTCC"
params.light_chain_motif    = "aGCgACgGGaGTtCAcagGTATACATGTTGCTGTGGTTGTCTG"
params.igk_idx              = 336
params.bcr_count_thresh     = 5



log.info """\
G C Re - F L O W!
Matsen, Victora Labs
Fred Hutchinson CRC, Seattle WA
Rockefeller University, New York NY.
================================
"""

/*
 * Import modules
 */

include {
    TRIM_COMBINE_MATES;
    DEMULTIPLEX_PLATES;
    DEMULTIPLEX_WELLS;
    SPLIT_HK as SPLIT_HEAVY;
    SPLIT_HK as SPLIT_LIGHT;
    COLLAPSE_RANK_PRUNE;
    MERGE_BCRS;
    PARTIS_ANNOTATION;
    PARTIS_WRANGLE;
    GCTREE;
    MERGE_RESULTS;
    NDS_LB_ANALYSIS;
    FITNESS_REGRESSION_ANALYSIS;
    MUTATIONS_ANALYSIS;
    INTERACTIVE_HEATMAPS;
  } from './modules.nf'


workflow BCR_COUNTS {

  take:
    filepair

  main:

    TRIM_COMBINE_MATES(filepair) | set { trimmed_ch }
    DEMULTIPLEX_PLATES(trimmed_ch, file(params.plate_barcodes)) \
      | transpose() | filter{ file(it[1]).size()>0 } \
      | set { dmplxd_plates_ch }

    DEMULTIPLEX_WELLS(dmplxd_plates_ch, file(params.well_barcodes)) \
      | transpose() | filter{ file(it[1]).size()>0 } \
      | set { dmplxd_wells_ch }

    SPLIT_HEAVY(
      dmplxd_wells_ch,
      params.heavy_chain_motif, "H"
    ) | filter{ file(it[1]).size()>0 } | set { heavy_ch }

    SPLIT_LIGHT(
      dmplxd_wells_ch,
      params.light_chain_motif, "K"
    ) | filter{ file(it[1]).size()>0 } | set { light_ch }

    heavy_ch.concat(light_ch) | COLLAPSE_RANK_PRUNE \
      | groupTuple(by:[0], sort:true) | MERGE_BCRS

  emit:
    MERGE_BCRS.out

}

workflow {

  Channel.fromPath(params.ngs_manifest)
    .splitCsv(header:true)
    .map{ row ->
      tuple(
        "$row.ngs_id",
        "$row.date",
        file("${params.reads_prefix}/${row.read1}"),
        file("${params.reads_prefix}/${row.read2}"),
      )
    } | BCR_COUNTS

  PARTIS_ANNOTATION(
    BCR_COUNTS.out, 
    file(params.partis_anno_dir)
  ) | set{ partis_anno_ch }

  PARTIS_WRANGLE(
    partis_anno_ch, 
    file(params.gc_metadata)
  ) | flatten() | set{ partis_wrangle_ch }

  GCTREE(
    partis_wrangle_ch, 
    file(params.hdag_sub), 
    file(params.hdag_mut), 
    file(params.dms_vscores),
    file(params.dms_sites),
    file("$projectDir/bin/gctree-tools.py"),
    file("$projectDir/bin/trees.py")
  ) | collect | set{gctree_ch}

  MERGE_RESULTS(gctree_ch)

  ranking_coeff_strategy_ch = Channel.of(
    "default", "naive_reversions_first", "naive_reversions_no_bp"
  )

  gctree_ch
    .map{it -> [it]}
    .combine(ranking_coeff_strategy_ch)
    .set{ gctree_rank_ch }

  NDS_LB_ANALYSIS(
    file("${projectDir}/notebooks/NDS-LB.ipynb"),
    file("${projectDir}/notebooks/utils/"),
    file(params.gc_metadata),
    gctree_rank_ch.combine(Channel.of(5, 20)) // svg scale
  )

  FITNESS_REGRESSION_ANALYSIS(
    file("${projectDir}/notebooks/fitness-regression.ipynb"),
    file("${projectDir}/notebooks/utils/"),
    file(params.gc_metadata),
    gctree_rank_ch
  )

  MUTATIONS_ANALYSIS(
    file("${projectDir}/notebooks/mutations.ipynb"),
    file("${projectDir}/notebooks/utils/"),
    file(params.gc_metadata),
    file(params.dms_vscores),
    file(params.dms_sites),
    file(params.chigy_hc_mut_rates),
    file(params.chigy_lc_mut_rates),
    file(params.pdb),
    gctree_rank_ch
  )

  INTERACTIVE_HEATMAPS(
    file("${projectDir}/notebooks/interactive-heatmaps.ipynb"),
    MUTATIONS_ANALYSIS.out | map{it -> tuple(it[2], it[3])}
  )

}
