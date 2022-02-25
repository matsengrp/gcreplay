# GC-Replay-NF

A Nextflow pipeline for running the analysis of germinal center replay 
(gcreplay) experimental analysis.


## Quickstart 

Install `Nextflow` by using the following command: 

    $ curl -s https://get.nextflow.io | bash 
    
Download the `Docker` Desktop, there exists several distibutions packaged for
various linux flavors

    $ curl -fsSL https://get.docker.com -o get-docker.sh && sudo sh get-docker.sh

Launch the pipeline execution with the following command: 

    $ nextflow run matsengrp/gcreplay -profile docker

Note: the [Dockerfile](docker/Dockerfile) contains all the required dependencies. 
Add the `-profile docker` to enable the containerised execution to the 
example command line shown below. 

**WHILE PRIVATE**

    $ git clone git@github.com:matsengrp/gcreplay.git && cd gcreplay
    $ git checkout 1_nextflow_infrastructure
    $ nextflow run nextflow/main.nf -profile docker -resume

```
N E X T F L O W  ~  version 21.04.3
Launching `main.nf` [hungry_mendel] - revision: fd014de395
G C Re - F L O W!
Matsen, Victora Labs
Fred Hutchinson CRC, Seattle WA
Rockefeller University, New York NY.
================================

executor >  slurm (94)
[d8/33858b] process > BCR_COUNTS:TRIM_COMBINE_MATES (5)      [100%] 8 of 8, cached: 8 ✔
[ba/75b8b9] process > BCR_COUNTS:DEMULTIPLEX_PLATES (2)      [100%] 8 of 8, cached: 8 ✔
[b9/7a4d21] process > BCR_COUNTS:DEMULTIPLEX_WELLS (105)     [100%] 109 of 109, cached: 109 ✔
[66/6441e4] process > BCR_COUNTS:SPLIT_HEAVY (9398)          [100%] 9444 of 9444, cached: 9444 ✔
[57/7d2930] process > BCR_COUNTS:SPLIT_LIGHT (9423)          [100%] 9444 of 9444, cached: 9444 ✔
[fb/a97f83] process > BCR_COUNTS:COLLAPSE_RANK_PRUNE (18799) [100%] 18888 of 18888, cached: 18888 ✔
[6e/86943b] process > BCR_COUNTS:MERGE_BCRS (2)              [100%] 8 of 8, cached: 8 ✔
[1d/f8a5af] process > PARTIS_ANNOTATION (8)                  [100%] 8 of 8, cached: 8 ✔
[d6/1c7de2] process > PARTIS_WRANGLE (7)                     [100%] 8 of 8 ✔
[b5/2b4b19] process > GCTREE (78)                            [100%] 86 of 86 ✔
Completed at: 24-Feb-2022 13:07:30
Duration    : 1h 23m 32s
CPU hours   : 247.1 (51.4% cached)
Succeeded   : 94
Cached      : 37'917
```

Above we see a description of the pipeline run. Total there were ~ 37K jobs, taking
sum total of ~ 250 CPU hours in the duration of about an hour and a half when run on
the Fred Hutchinson HPC. 

## Pipeline Description

(Currently) This pipeline processes reads from a paired-end MiSeq run and generally 
performs the following steps; 
(1) trim and combine paired end files 
(2) demultiplex both plates _and_ wells
(3) Split Heavy and light chain reads per well
(4) Collapse identical sequences, while retaining and ordering rank of each sequence per well
(5) Prune to keep the top N sequences observed from each well
(6) Merge the results, formatting for [partis](TODO) annotation.
(7) [partis](TODO) annotation
(8) curation, cleaning, and merging of Heavy and Light Chains
(9) [gctree](TODO) lineage inference using HDAG.


## Input files

The CalliNGS-NF pipeline needs as the input following files:
* Paired-end MiSeq reads, `*.fastq` (currently must be decompressed)
* A barcode file for demultiplexing plates (txt)
* A barcode file for demultiplexing individual plate wells from FACS single sort
* Germline sequence information for partis annotation

## Pipeline parameters

#### `--manifest` 
   
A CSV describing 
1. "sample_id" A unique identifier for and set of paired miseq NGS files (a "PR")
2. "key_file" A Key File describing the Germinal Centers expected in each well
2. "date" Date the experiment was sequenced 
3. "read1" & "read2" Paths to each of the paried end fastq files - relative to the project directory (where the main.nf file is located)

#### `--plate_barcodes`

* Specifies the barcodes (at the **beginning** of read) for demultiplexing
individual plates from single pair of MiSeq Raw NGS Fastq files

* formatted with two, space-deliminated, columns. The first column
should be the plate ID, and the second the oligonucleotide barcode  

#### `--well_barcodes`

* Specifies the barcodes (at the **end** of read) for demultiplexing
individual wells from each of the plate-demultiplexed NGS Fastq files.

* formatted with two, space-deliminated, columns. The first column
should be the well ID, and the second the oligonucleotide barcode  

#### `--partis_anno_dir`

* Path to the germline sequences needed for partis annotation.
- relative to the project directory (where the main.nf file is located)

#### `--top_n_rank`

* This parameter is used after


#### `results`

* where you would like the final results to be places

* a relative path (does not need to exist)

    
## Pipeline results

TODO

## Schematic Outline

<img src="./data/images/dag.svg">

## Requirements 

TODO

required software components reported in the following section. See the included 
[Dockerfile](docker/Dockerfile) for the configuration details.
 
## Components 

TODO

