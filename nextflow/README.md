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

## Pipeline Description

(Currently) This pipeline processes reads from a paired-end MiSeq run and generally 
performs the following steps; 
(1) trim and combine paired end files 
(2) demultiplex both plates _and_ wells
(3) Split Heavy and light chain reads per well
(4) Collapse identical sequences, while retaining and ordering rank of each sequence per well
(5) Prune to keep the top N sequences observed from each well
(6) Merge the results, formatting for [partis]() annotation.


## Input files

The CalliNGS-NF pipeline needs as the input following files:
* Paired-end MiSeq reads, `*.fastq` (currently must be decompressed)
* A barcode file for demultiplexing plates (txt)
* A barcode file for demultiplecing

## Pipeline parameters

#### `--reads` 
   
* Specifies the location of the reads FASTQ file(s).

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
