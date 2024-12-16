/*
 * This Source Code Form is subject to the terms of the GNU GENERAL PUBLIC LICENCE
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
params.reads_prefix     = "$projectDir"
params.manifest         = "data/test/manifest.csv"
// params.metadata         = "data/metadata/metadata.csv"
params.plate_barcodes   = "data/barcodes/plateBC.txt"
params.well_barcodes    = "data/barcodes/96FBC.txt"
params.partis_anno_dir  = "$projectDir/data/partis_annotation/germlines"
params.results          = "$projectDir/results/"
params.hdag_sub         = "data/mutability/MK_RS5NF_substitution.csv"
params.hdag_mut         = "data/mutability/MK_RS5NF_mutability.csv"
params.dms_vscores      = "data/dms/final_variant_scores.csv"
params.dms_sites        = "data/dms/CGGnaive_sites.csv"
params.igk_idx          = 336
params.bcr_count_thresh = 5



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
  } from './modules.nf'


workflow BCR_COUNTS {

  take:
    filepair

  main:

    TRIM_COMBINE_MATES(filepair) | set { trimmed_ch }
    DEMULTIPLEX_PLATES(trimmed_ch, file("$params.plate_barcodes")) \
      | transpose() | filter{ file(it[2]).size()>0 } \
      | set { dmplxd_plates_ch }

    DEMULTIPLEX_WELLS(dmplxd_plates_ch, file("$params.well_barcodes")) \
      | transpose() | filter{ file(it[2]).size()>0 } \
      | set { dmplxd_wells_ch }

    SPLIT_HEAVY(
      dmplxd_wells_ch,
      "aGCgACgGGaGTtCAcagACTGCAACCGGTGTACATTCC", "H"
    ) | filter{ file(it[2]).size()>0 } | set { heavy_ch }

    SPLIT_LIGHT(
      dmplxd_wells_ch,
      "aGCgACgGGaGTtCAcagGTATACATGTTGCTGTGGTTGTCTG", "K"
    ) | filter{ file(it[2]).size()>0 } | set { light_ch }

    heavy_ch.concat(light_ch) | COLLAPSE_RANK_PRUNE \
      | groupTuple(by:[0,1], sort:true) | MERGE_BCRS

  emit:
    MERGE_BCRS.out

}

workflow {

  Channel.fromPath(params.manifest)
    .splitCsv(header:true)
    .map{ row ->
      tuple(
        "$row.sample_id",
        file("${params.reads_prefix}/${row.key_file}"),
        "$row.date",
        file("${params.reads_prefix}/${row.read1}"),
        file("${params.reads_prefix}/${row.read2}"),
      )
    } | BCR_COUNTS

  // PARTIS_ANNOTATION(BCR_COUNTS.out) \
  //   | PARTIS_WRANGLE | flatten() | set{partis_wrangle_ch}

  // GCTREE(
  //   partis_wrangle_ch, 
  //   file("$params.hdag_sub"), 
  //   file("$params.hdag_mut"), 
  //   file("$params.dms_vscores"),
  //   file("$params.dms_sites")
  // ) | collect | set{gctree_ch}
  
  // MERGE_RESULTS(gctree_ch)
  
  // // Channel for metadata file as value
  // metadata_ch = Channel.value(file("${projectDir}/${params.metadata}"))

  // // Channel for ranking coefficients
  // ranking_coeff_ch = Channel.of("default", "naive_reversions_first", "naive_reversions_no_bp")

  // // [[PR_dir, PR_dir2], metadata, ranking_coeff]
  // gctree_meta_rank_ch = gctree_ch
  //   .map{it -> [it]}
  //   .combine(metadata_ch)
  //   .combine(ranking_coeff_ch)
  
  // // run iterations of scale, and rank through the NDS_LB_ANALYSIS process
  // gctree_meta_rank_ch.combine(Channel.of(5, 20)) | NDS_LB_ANALYSIS

}
