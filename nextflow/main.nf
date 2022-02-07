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
 * TODO finish all authors
 */

/* 
 * Enable DSL 2 syntax
 */
nextflow.enable.dsl = 2

/*
 * Define the default parameters - example data get's run by default
 */ 

params.manifest = "$baseDir/data/test/manifest.csv"
if (params.manifest != "$baseDir/data/test/manifest.csv")
    params.reads_prefix = "$launchDir"
else
    params.reads_prefix = "$baseDir"

params.plate_barcodes   = "$baseDir/data/barcodes/test_plateBC.txt"
params.well_barcodes    = "$baseDir/data/barcodes/test_96FBC.txt"

// if we're not using the default test, make the filespaths in the
// manifest relative to the launch directory (assuming files are now local)

params.partis_anno_dir  = "$baseDir/data/partis_annotation/"
params.results          = "$launchDir/results/"


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
  } from './modules.nf' 


workflow BCR_COUNTS {

  take: 
    filepair

  main:

    TRIM_COMBINE_MATES(filepair)
    DEMULTIPLEX_PLATES(TRIM_COMBINE_MATES.out) \
      | transpose() | filter{ file(it[1]).size()>0 } \
      | set { dmplxd_plates_ch }
    
    DEMULTIPLEX_WELLS(dmplxd_plates_ch) \
      | transpose() | filter{ file(it[1]).size()>0 } \
      | set { dmplxd_wells_ch }

    SPLIT_HEAVY( 
      dmplxd_wells_ch, 
      "aGCgACgGGaGTtCAcagACTGCAACCGGTGTACATTCC", "H"  
    )

    SPLIT_LIGHT( 
      dmplxd_wells_ch, 
      "aGCgACgGGaGTtCAcagGTATACATGTTGCTGTGGTTGTCTG", "K"  
    )

    SPLIT_HEAVY.out.mix(SPLIT_LIGHT.out) | COLLAPSE_RANK_PRUNE \
        | groupTuple() | MERGE_BCRS

  emit:
    MERGE_BCRS.out

}

workflow {

  /*
   * WHAT I WOULD LIKE TO BE DOING
   * I would like to fire off a BCR_COUNTS workflow
   * for each file in the manifest
   */

  // TODO add date to sequence id
  Channel.fromPath(params.manifest)
    .splitCsv(header:true)
    .map{ row -> 
      tuple(
        "$row.sample_id",
        "$row.date",
        file("${params.reads_prefix}/${row.read1}"),
        file("${params.reads_prefix}/${row.read2}"),
      )
    } | BCR_COUNTS

  // Step 2
  PARTIS_ANNOTATION(BCR_COUNTS.out)

  // Step 3
  // | PREP_ANNOTATION

  // Step 4
  // | GCTREE

  // Step 5
  // | DATABASE_WRANGLE

}






