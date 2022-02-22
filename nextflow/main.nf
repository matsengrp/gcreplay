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

// csv with columns 
params.manifest = "$baseDir/data/test/manifest.csv"

// if we're not using the default test, make the filespaths in the
// manifest relative to the launch directory (assuming files are now local)
if (params.manifest != "$baseDir/data/test/manifest.csv")
    params.reads_prefix = "$launchDir"
else
    params.reads_prefix = "$baseDir"

// all plate barcodes.
params.plate_barcodes   = "$baseDir/data/barcodes/test_plateBC.txt"

// all the 96 well barcodes. 
params.well_barcodes    = "$baseDir/data/barcodes/test_96FBC.txt"


// temp - where do you want the partis annotation, specfically.
params.partis_anno_dir  = "$baseDir/data/partis_annotation/"

// the directory you would like to plave all the results
params.results          = "$launchDir/results/"

// keep n lines of a sequence collapsed fasta.
// either this or the parameter below must be zero
// TODO do the multiplication by 2 yourself in the process script block
params.top_n_rank       = 6

// keep all sequences above
// TODO implement 
params.n_threshhold     = 0


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
  } from './modules.nf' 


workflow BCR_COUNTS {

  take: 
    filepair

  main:

    TRIM_COMBINE_MATES(filepair)
    DEMULTIPLEX_PLATES(TRIM_COMBINE_MATES.out) \
      | transpose() | filter{ file(it[2]).size()>0 } \
      | set { dmplxd_plates_ch }
    
    DEMULTIPLEX_WELLS(dmplxd_plates_ch) \
      | transpose() | filter{ file(it[2]).size()>0 } \
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
        | groupTuple(by:[0,1]) | MERGE_BCRS

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
        file("${params.reads_prefix}/${row.key_file}"),
        "$row.date",
        file("${params.reads_prefix}/${row.read1}"),
        file("${params.reads_prefix}/${row.read2}"),
      )
    } | BCR_COUNTS

  // Step 2
  PARTIS_ANNOTATION(BCR_COUNTS.out) | PARTIS_WRANGLE | view()
  

  // Step 3
  // | PREP_ANNOTATION

  // Step 4
  // | GCTREE

  // Step 5
  // | DATABASE_WRANGLE

}






