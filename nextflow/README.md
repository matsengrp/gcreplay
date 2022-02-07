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
Launching `main.nf` [deadly_legentil] - revision: 2ed283eecb
G C Re - F L O W!
Matsen, Victora Labs
Fred Hutchinson CRC, Seattle WA
Rockefeller University, New York NY.
================================

executor >  slurm (37871)
[fa/5eb808] process > BCR_COUNTS:TRIM_COMBINE_MATES (6)      [100%] 8 of 8, cached: 3 ✔
[b4/4d007b] process > BCR_COUNTS:DEMULTIPLEX_PLATES (8)      [100%] 8 of 8, cached: 3 ✔
[24/35cb6a] process > BCR_COUNTS:DEMULTIPLEX_WELLS (107)     [100%] 109 of 109, cached: 32 ✔
[9a/6d6604] process > BCR_COUNTS:SPLIT_HEAVY (9414)          [100%] 9442 of 9442 ✔
[dc/e940b3] process > BCR_COUNTS:SPLIT_LIGHT (9429)          [100%] 9442 of 9442 ✔
[46/195631] process > BCR_COUNTS:COLLAPSE_RANK_PRUNE (18856) [100%] 18884 of 18884 ✔
[36/197533] process > BCR_COUNTS:MERGE_BCRS (6)              [100%] 8 of 8 ✔
[f3/9ea99d] process > PARTIS_ANNOTATION (8)                  [100%] 8 of 8 ✔
Completed at: 06-Feb-2022 22:21:22
Duration    : 44m 27s
CPU hours   : 24.1 (2% cached)
Succeeded   : 37'871
Cached      : 38
```

## Pipeline Description

(Currently) This pipeline processes reads from a paired-end MiSeq run and generally 
performs the following steps; 
(1) trim and combine paired end files 
(2) demultiplex both plates _and_ wells
(3) Split Heavy and light chain reads per well
(4) Collapse identical sequences, while retaining and ordering rank of each sequence per well
(5) Prune to keep the top N sequences observed from each well
(6) Merge the results, formatting for [partis]() annotation.
(7) partis annotation


## Input files

The CalliNGS-NF pipeline needs as the input following files:
* Paired-end MiSeq reads, `*.fastq` (currently must be decompressed)
* A barcode file for demultiplexing plates (txt)
* A barcode file for demultiplecing

## Pipeline parameters

#### `--manifest` 
   
* TODO

#### `plate_barcodes`

* Specifies the barcodes (at the **beginning** of read) for demultiplexing
individual plates from single pair of MiSeq Raw NGS Fastq files

* formatted with two, space-deliminated, columns. The first column
should be the plate ID, and the second the oligonucleotide barcode  

#### `well_barcodes`

* Specifies the barcodes (at the **end** of read) for demultiplexing
individual wells from each of the plate-demultiplexed NGS Fastq files.

* formatted with two, space-deliminated, columns. The first column
should be the well ID, and the second the oligonucleotide barcode  

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
