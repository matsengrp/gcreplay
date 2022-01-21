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

params.reads        = "$baseDir/data/test/*_r{1,2}.fastq"
params.plate_bc     = "$baseDir/data/test/plateBC.txt"
params.well_bc      = "$baseDir/data/test/96FBC.txt"
params.results      = "results"


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

workflow {
  reads_ch = Channel.fromFilePairs(params.reads)
  plate_bc = Channel.fromPath(params.plate_bc)
  well_bc = Channel.fromPath(params.well_bc)

  TRIM_COMBINE_MATES(reads_ch)
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
  | collect | MERGE_BCRS 
}








