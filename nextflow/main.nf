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
// if we're not using the default test, make the filespaths in the
// manifest relative to the launch directory (assuming files are now local)
if (params.manifest != "$baseDir/data/test/manifest.csv")
    params.reads_prefix = "$launchDir"
else // otherwise we need them relative to remote repo (baseDir)
    params.reads_prefix = "$baseDir"

params.plate_bc     = "$baseDir/data/barcodes/plateBC.txt"
params.well_bc      = "$baseDir/data/barcodes/96FBC.txt"
params.results      = "$launchDir/results/"


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
  } from './modules.nf' 


workflow BCR_COUNTS {

  take: 
    filepair

  main:

    plate_bc = Channel.fromPath(params.plate_bc).first()
    well_bc = Channel.fromPath(params.well_bc).first()

    TRIM_COMBINE_MATES(filepair)
    DEMULTIPLEX_PLATES(TRIM_COMBINE_MATES.out, plate_bc) \
      | flatten() | filter{ file(it).size()>0 } \
      | combine(well_bc) | set { dmplxd_plates_ch }
    DEMULTIPLEX_WELLS(dmplxd_plates_ch) \
      | flatten() | filter{ file(it).size()>0 } \
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
    | collect | set { all_ranked_ch }
    MERGE_BCRS(filepair, all_ranked_ch)

  emit:
    MERGE_BCRS.out

}

workflow {

  Channel.fromPath(params.manifest)
    .splitCsv(header:true)
    .map{ row -> 
      tuple(
        "$row.sample_id",
        file("${params.reads_prefix}/${row.read1}"),
        file("${params.reads_prefix}/${row.read2}"),
      )
    } | BCR_COUNTS
    
  BCR_COUNTS.out | view()

  // now for annotation

}






