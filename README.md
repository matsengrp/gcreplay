# GC replay

This repository contains the code and data for the analysis of B cell receptor (BCR) sequences from germinal center (GC) B cells in mice.

The repository is currently broken into three main directories: [nextflow](nextflow/) contains the nextflow pipeline which processes the raw data and generates the BCR trees, [analysis](analysis/) contains a set of jupyter notebooks used in the downstream analysis of these counts, and finally, [passenger](passenger/) contains a set of scripts used to analyze the passenger mouse data. See the READMEs in each directory for more information on each.


## Description of metadata file metadata.csv

The metadata file `metadata.csv` contains information about the samples used in the analysis. The columns are as follows:

- "uid": (_unique identifier_) This is in the format of D<imm_duration>_M<mouse>_GC<gc>. 
- "ngs_id": (_Parallel Replay_) reference to associated sequencing run specified in [ngs_manifest.csv](ngs_manifest.csv).
- imm_duration: (_immunization duration_) This is the number of weeks after immunization that the sample was taken.
- mouse: (_mouse number_) This is the unique identifier for the mouse.
- gc: (_germinal center number_) This is the unique identifier for the germinal center.
- strain: (_mouse strain_) This is the strain of the mouse.
- node: (_lymph node number_) This identifies the specific lymph node location in the mouse.
- cell_type: (_cell type_) This is the type of cell that was sequenced.
- plate: (_plate barcode_) This is the barcode of the plate that the sample was sequenced on.
- hc_barcode : (_heavy chain barcode_) This is the barcode of the heavy chain.
- lc_barcode : (_light chain barcode_) This is the barcode of the light chain.
- row : (_well rows_) Which rows of the 96-well plate contain the sample.
- col : (_well columns_) Which columns of the 96-well plate contain the sample.

In september of 2024, and December 2024, the mouse and GC identifiers were changed to be more informative. The old identifiers are stored in [reindex_mapping.csv](reindex_mapping.csv) for reference.


### other notes on non-GC samples.

- the w10 (week 10) samples used lymph node (LN) sampling rather than GC sampling because there weren't any GCs big enough to find (GC number for these is arbitrary, at the moment it's just sampling order)                                                                                      

- LMP2a mice had LMP2A protein from virus inserted into one of the D regions of the VH locus
    - So the B cell expresses LMP2A and chigy-H from the other VH locus
    - Effect is that LMP2A providing weak constitutively active B cell signaling (tonic signal)
    - We thought that this might provide B cells that have nonsense mutations with enough B cell signals to survive
    - If this was the case it would’ve been consistent with a hypothesis that there might be a checkpoint from leaving the dark zone to enter the right zone (ie having receptors on surface)
    - but we found no evidence of that, meaning LMP2A GCs didn’t increase B cells with stop codon
    - In fact it might have increased B cells with affinity increases (more efficient affinity maturation?) but don’t remember if there was statistical significance# GC-Replay-NF

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


## Pipeline results

The primary results produced by the main pipeline will be stored under the specified `--results`
path. Beneath this you'll find a few different subdirectories containing the intermediate files
from individual steps of the pipeline. The primary outputs consist of the individual inferred trees
for each clonal family. The most relevant outputs that we track in this repository are:


### [results/merged-results](results/archive/2024-10-09-rank-coeff-renumbered/merged-results)

Here, you'll find the summary of all tree nodes across all trees, and all observed BCR's in the form of two tabular csv files. 

[observed-seqs.csv](results/merged-results/observed-seqs.csv) contains the observed BCR sequences, with the following columns:

- ID_HK: unique identifier for the BCR
- well: well from which the BCR was identified, in the format of <row><column>
- HK_key_plate: barcode of the plate
- HK_key_mouse: mouse identifier (according to the metadata "mouse" col)
- HK_key_gc: germinal center identifier (according to the metadata "gc" col)
- HK_key_node: lymph node identifier (according to the metadata "node" col)
- HK_key_cell_type: cell type identifier (according to the metadata "cell_type" col)
- aa_substitutions_IMGT: amino acid substitutions in IMGT format

The next three columns are the inferred delta values for each of the three metrics, according to the main pipeline params.dms_vscores, TODO point to Kd pipeline?

- delta_bind: delta binding energy, as inferred through additive Kd model
- delta_expr: delta expression, as inferred through additive expression values of individual mutations
- delta_psr: delta polyspecificity, as inferred through additive polyspecificity values of individual mutations

The next four columns are the number of mutations in the heavy and light chains, and the number of amino acid mutations in the heavy and light chains.

- n_nt_mutations_HC: number of nucleotide mutations in the heavy chain
- n_nt_mutations_LC: number of nucleotide mutations in the light chain
- IgH_nt_mutations: nucleotide mutations in the heavy chain, formatted like `<wt><sequential-site>-<mut>` 
- IgK_nt_mutations: nucleotide mutations in the light chain, formatted like `<wt><sequential-site>-<mut>`
- n_aa_mutations_HC: number of amino acid mutations in the heavy chain
- n_aa_mutations_LC: number of amino acid mutations in the light chain
- IgH_aa_mutations: Same as 'aa_substitutions_IMGT', but just the IgH mutations
- IgK_aa_mutations: Same as 'aa_substitutions_IMGT', but just the IgK mutations

The next two columns are the isotypes of the heavy and light chains. Isotypes are inferred via fuzzy motif matching

- isotype_HC: isotype of the heavy chain
- isotype_LC: isotype of the light chain

The next two columns are the unique identifiers for the heavy and light chains - they are essentially the same as the ID_HK, but one for each chain

- ID_HC: unique identifier for the heavy chain
- ID_LC: unique identifier for the light chain

The following columns are partis annotations: TODO ask duncan to annotate?

- Productive_HC: 
- Productive_LC:
- V_HC:
- V_LC:
- D_HC:
- D_LC:
- J_HC:
- J_LC:
- AAjunction_HC:
- AAjunction_LC:
- locus_HC:
- locus_LC:

Other misc annotations regarding the barcode and plate locations the cell derived from

- miseq_plate_HC:
- miseq_plate_LC:
- barcode_HC:
- barcode_LC:
- row_HC: The row of the 96 well plate from which this cell was derived
- row_LC: The row of the 96 well plate from which this cell was derived
- column_HC: The column of the 96 well plate from which this cell was derived
- column_LC: The column of the 96 well plate from which this cell was derived
- chain_HC: ?
- chain_LC: ?
- rank_HC: The rank of the cell within the well. ranked from 1 - N seqs that met the bcr count threshold
- rank_LC: The rank of the cell within the well. ranked from 1 - N seqs that met the bcr count threshold
- counts_HC: The number of times this cell was observed in the well
- counts_LC: The number of times this cell was observed in the well
- seq_nt_length_HC: The length of the nucleotide sequence of the heavy chain
- seq_nt_length_LC: The length of the nucleotide sequence of the light chain
- seq_aa_length_HC: The length of the amino acid sequence of the heavy chain
- seq_aa_length_LC: The length of the amino acid sequence of the light chain
- fasta_header_HC: The fasta header of the heavy chain
- fasta_header_LC: The fasta header of the light chain
- fasta_seq_HC: The fasta sequence of the heavy chain
- fasta_seq_LC: The fasta sequence of the light chain
- partis_sequence_HC: The partis sequence of the heavy chain
- partis_sequence_LC: The partis sequence of the light chain
- seq_aa_HC: The amino acid sequence of the heavy chain
- seq_aa_LC: The amino acid sequence of the light chain
- seq_nt_HC: The nucleotide sequence of the heavy chain
- seq_nt_LC: The nucleotide sequence of the light chain


### [results/gctrees](results/archive/2024-10-09-rank-coeff-renumbered/gctrees), 

Here, you'll find subdirectories containing all intermediate input and output files from running `gctree`. Each directory is named in the format `PR<PR>-<mouse>-<node>-<gc>`
where each of the `<fields>` can be matched to the [metadata](../metadata.csv) file to identify the mouse. TODO: this should be updated to just be the <uid>.

Each gctree subdirectory now contains results for each three ranking coeff approaches used in the manuscript. For example, in the [results/gctrees/PR1.01-20-LB-55-GC](results/gctrees/PR1.01-20-LB-55-GC) directory, after execution of the pipeline, you'll find the following files:
that were used as input to gctree:
```
├── abundances.csv
├── annotated-PR-1-01-20-LB-55-GC.csv
├── annotated-PR-1-01-20-LB-55-GC.fasta
├── annotated-PR-1-01-20-LB-55-GC.idmap
├── annotated-PR-1-01-20-LB-55-GC.isotypemap
├── default/
├── naive_reversions_first/
├── naive_reversions_no_bp/
├── observed_seqs.csv
└── outfile
```
where each of the "default", "naive_reversions_first", "naive_reversions_no_bp" strategy directories looks like

```
nextflow/results/2024-09-30-rank-coeff-renumbered/gctrees/PR1.01-20-LB-55-GC/default
├── delta_bind.fasta
├── delta_bind.svg
├── delta_expr.fasta
├── delta_expr.svg
├── delta_psr.fasta
├── delta_psr.svg
├── gctree_default.inference.log
├── gctree.forest_summary.log
├── gctree.inference.1.nk
├── gctree.inference.1.p
├── gctree.inference.1.svg
├── gctree.inference.abundance_rank.svg
├── gctree.inference.parsimony_forest.p
├── gctree.p
├── LBI.fasta
├── LBI.svg
├── LBR.fasta
├── LBR.svg
└── node_data.csv
```

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

* formatted with two, space-delimited, columns. The first column
should be the well ID, and the second the oligonucleotide barcode  

#### `--partis_anno_dir`

* Path to the germline sequences needed for partis annotation.
- relative to the project directory (where the main.nf file is located)

#### `--top_n_rank`

* This parameter is used after


#### `--results`

* where you would like the final results to be places

* a relative path (does not need to exist)




## Schematic Outline

<img src="./data/images/dag.svg">

## Requirements 

required software components reported in the following section. See the included 
[Dockerfile](docker/Dockerfile) for the configuration details.

# Beast pipeline

[beast.nf](beast.nf) runs [BEAST (v1)](https://beast.community/) on a set of naive and observed BCR sequences 
from a single germinal center to infer time trees for each clonal family.

In more detail, the pipeline prepares the xml files from [beastgen](https://beast.community/beastgen) 
and a specified template 
(pre-configured templates can be found in the 
[data/beast/beast_templates](data/beast/beast_templates) directory). 
By default, the pipeline uses the [skyline histlog](/data/beast/beast_templates/skyline_histlog.template.patch) template.
The pipeline then patches the xml to fix the naive sequence in time using the [beast_template_root_fix.py script](bin/beast_template_root_fix.py)
before running BEAST on the patched xml files.


## Quick start

There's a testing set of sequences in the [data/beast/](data/beast/) directory. To run the beast pipeline on this data (with a small number of MCMC iterations for quick execution), you can use the following command:

```
nextflow run beast/main.nf --chain_length 1000 --log_every 100 -profile docker -resume
```

## Pipeline parameters

* `--seqs` is a string parameter specifying a file path pattern to the fasta files containing the naive and observed BCR sequences. The beast pipeline will be run on each of the files matching this pattern.

* `--beast_template` is a string parameter specifying the path to the beast template file. The pipeline will use this template to generate the xml files for the beast runs.

* `--results` is a string parameter specifying the path to the directory where the results of the beast runs will be stored.

* `--chain_length` is an integer parameter specifying the number of MCMC iterations to run the beast for.

* `--log_every` is an integer parameter specifying the interval at which MCMC-step trees are recorded.

* `--convert_to_ete` is a boolean parameter specifying whether to convert the beast trees to ete trees.

* `--dms_vscores` is the url to the dms variant scores for adding phenotypes to the ete converted trees.

* `--dms_sites` is the url to the dms sites for adding phenotypes to the ete converted trees.

* `--burn_frac` is the fraction of the chain to discard as burn-in instead of converting to ete.

* `--save_pkl_trees` is a boolean parameter specifying whether to save the ete trees as pickle files. This can be very memory intensive when there are many logged tree iterations for each tree.

