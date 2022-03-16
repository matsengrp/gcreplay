# GC-Replay-NF

A Nextflow pipeline for running the analysis of germinal center replay 
(gcreplay) experimental analysis.


## Quickstart

Install `Nextflow` by using the following command:

    $ curl -s https://get.nextflow.io | bash

Download the `Docker` Desktop, there exists several distibutions packaged for
various linux flavors

    $ curl -fsSL https://get.docker.com -o get-docker.sh && sudo sh get-docker.sh

Note: the [Dockerfile](docker/Dockerfile) contains all the required dependencies.
Add the `-profile docker` to enable the containerised execution to the
example command line shown below.

Launch the pipeline execution with the following command:

    $ git clone git@github.com:matsengrp/gcreplay.git && cd gcreplay
    $ nextflow run nextflow/main.nf -profile docker -resume

## Obtain results

if you haven't already cloned, simply [install git lfs](https://git-lfs.github.com/)

    $ sudo apt install git-lfs
    $ git clone git@github.com:matsengrp/gcreplay.git && cd gcreplay

If you've already cloned, install git lfs, move into the repo and

    $ git lfs install

The final observed BCR and gctree results will then be in the [results](nextflow/results/) directory,
named by date of the pipeline execution producing each respective set
of results.

## Rhino Cluster

If you have obtained the necessary NGS files, and would like to run the full pipeline
on the fred hutch rhino cluster, you might run a script like the following

```
#!/bin/bash

set -e
source /app/lmod/lmod/init/profile

module load nextflow
module load Singularity
export PATH=$SINGULARITYROOT/bin/:$PATH

nextflow run main.nf \
        --manifest "data/_ignore/miseq/manifest.csv" \
        --results "results/$(date -I)" \
        -work-dir /fh/scratch/delete30/matsen_e/jgallowa/work/ \
        -profile fred_hutch_rhino \
        -resume
```

Running the full pipeline will
print a summary to stdout


```
N E X T F L O W  ~  version 21.04.3

G C Re - F L O W!
Matsen, Victora Labs
Fred Hutchinson CRC, Seattle WA
Rockefeller University, New York NY.
================================

executor >  slurm (98)
[e7/4c53a6] process > BCR_COUNTS:TRIM_COMBINE_MATES (1)      [100%] 8 of 8, cached: 8 ✔
[53/be72b1] process > BCR_COUNTS:DEMULTIPLEX_PLATES (2)      [100%] 8 of 8, cached: 8 ✔
[37/0d8693] process > BCR_COUNTS:DEMULTIPLEX_WELLS (96)      [100%] 109 of 109, cached: 109 ✔
[0b/2cbb87] process > BCR_COUNTS:SPLIT_HEAVY (9348)          [100%] 9444 of 9444, cached: 9444 ✔
[f4/204e1c] process > BCR_COUNTS:SPLIT_LIGHT (9444)          [100%] 9444 of 9444, cached: 9444 ✔
[13/5cb9a6] process > BCR_COUNTS:COLLAPSE_RANK_PRUNE (17736) [100%] 17759 of 17759, cached: 17759 ✔
[90/612bbc] process > BCR_COUNTS:MERGE_BCRS (7)              [100%] 8 of 8, cached: 8 ✔
[1a/cda920] process > PARTIS_ANNOTATION (8)                  [100%] 8 of 8, cached: 8 ✔
[9f/79539f] process > PARTIS_WRANGLE (7)                     [100%] 8 of 8 ✔
[cf/56a87a] process > GCTREE (1)                             [100%] 89 of 89 ✔
[f1/d4f884] process > MERGE_RESULTS                          [100%] 1 of 1 ✔
Completed at: 15-Mar-2022 22:47:27
Duration    : 1h 30m 15s
CPU hours   : 443.7 (70.9% cached)
Succeeded   : 98
Cached      : 36'788
```

Above we see a description of the pipeline run. Total there were ~ 40K jobs, taking
sum total of ~ 400 CPU hours in the duration of about an hour and a half when run on
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
(10) Merge the results


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

