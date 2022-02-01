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

//params.reads            = "$baseDir/data/test/test_PR1.6_r{1,2}.fastq"
//params.plate_barcodes   = "$baseDir/data/barcodes/plateBC.txt"
//params.well_barcodes    = "$baseDir/data/barcodes/96FBC.txt"
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
    //tuple val(sample_id), path(read1), path(read2)

  main:

    TRIM_COMBINE_MATES(filepair)
    //TRIM_COMBINE_MATES([sample_id, read1, read2])
    DEMULTIPLEX_PLATES(TRIM_COMBINE_MATES.out) \
      | transpose() | filter{ file(it[1]).size()>0 } \
      | set { dmplxd_plates_ch }
      //| view()
    
    DEMULTIPLEX_WELLS(dmplxd_plates_ch) \
      | transpose() | filter{ file(it[1]).size()>0 } \
      | set { dmplxd_wells_ch }
    //dmplxd_wells_ch.view()

    SPLIT_HEAVY( 
      dmplxd_wells_ch, 
      "aGCgACgGGaGTtCAcagACTGCAACCGGTGTACATTCC", "H"  
    )

    SPLIT_LIGHT( 
      dmplxd_wells_ch, 
      "aGCgACgGGaGTtCAcagGTATACATGTTGCTGTGGTTGTCTG", "K"  
    )

    //SPLIT_HEAVY.out.mix(SPLIT_LIGHT.out) | COLLAPSE_RANK_PRUNE | view()
    
    SPLIT_HEAVY.out.mix(SPLIT_LIGHT.out) | COLLAPSE_RANK_PRUNE \
        | groupTuple() | MERGE_BCRS
    //| view()
    //| collect | set { all_ranked_ch }
    //MERGE_BCRS(all_ranked_ch)

  emit:
    MERGE_BCRS.out

}

workflow {

  /*
   * WHAT I WOULD LIKE TO BE DOING
   * I would like to fire off a BCR_COUNTS workflow
   * for each file in the manifest
   */

  Channel.fromPath(params.manifest)
    .splitCsv(header:true)
    .map{ row -> 
      tuple(
        "$row.sample_id",
        file("${params.reads_prefix}/${row.read1}"),
        file("${params.reads_prefix}/${row.read2}"),
      )
    } | BCR_COUNTS


  /*
   * FOR NOW
   * just call each one separately, I guess 
   * b/c I don't think calliong a workflow 
   * per record ina csv is possible, we may need
   * to just effectively group the keys, file
   * groupings throughout the BCR_COUNTS workflow
   */

  // Step 1
  //Channel.fromFilePairs(params.reads)
  //  .map { key, files ->
  //      tuple( 
  //          key,
  //          file(files[0]),
  //          file(files[1])
  //      )
  //  } | BCR_COUNTS


  // Step 2
  PARTIS_ANNOTATION(BCR_COUNTS.out)

  //BCR_COUNTS.out | view()
  // now for annotation

}






