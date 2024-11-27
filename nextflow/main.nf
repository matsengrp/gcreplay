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


// TODO let's move all these to the config file, huh?
// TODO we could just forget the example data, make them clone
// TODO for testing we should just make stub's
/*
 * Define the default parameters - example data get's run by default
 */

// rename to "manifest_path_prefix"
params.reads_prefix     = "$projectDir"
params.manifest         = "data/test/manifest.csv"
params.plate_barcodes   = "data/barcodes/plateBC.txt"
params.well_barcodes    = "data/barcodes/96FBC.txt"
params.partis_anno_dir  = "$projectDir/data/partis_annotation/germlines"
params.results          = "$projectDir/results/"
params.hdag_sub         = "data/mutability/MK_RS5NF_substitution.csv"
params.hdag_mut         = "data/mutability/MK_RS5NF_mutability.csv"
params.dms_vscores      = "data/dms/final_variant_scores.csv"
params.dms_sites        = "https://raw.githubusercontent.com/jbloomlab/Ab-CGGnaive_DMS/main/data/CGGnaive_sites.csv"
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

    heavy_ch.mix(light_ch) | COLLAPSE_RANK_PRUNE \
        | groupTuple(by:[0,1]) | MERGE_BCRS

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

  PARTIS_ANNOTATION(BCR_COUNTS.out) \
    | PARTIS_WRANGLE | flatten() | set{partis_wrangle_ch}

  GCTREE(partis_wrangle_ch, file("$params.hdag_sub"), file("$params.hdag_mut"), file("$params.dms_vscores")) \
    | collect | MERGE_RESULTS

}
