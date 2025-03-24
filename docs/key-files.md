# GCReplay Pipeline Results Documentation

This document provides an overview of the key files in the [gcreplay repository](https://github.com/matsengrp/gcreplay).

### [results/gc_metadata.csv](https://github.com/matsengrp/gcreplay/blob/main/gc_metadata.csv)

Contains metadata for each germinal center (GC) analyzed in the study. It includes information such as the day of analysis, mouse number, and GC identifier (uid). 

### [results/ngs_manifest.csv](https://github.com/matsengrp/gcreplay/blob/main/ngs_manifest.csv)
Serves as a manifest for the Next Generation Sequencing (NGS) data, detailing the sequencing runs and their corresponding metadata. It includes information such as the sequencing date, sample identifiers, and other relevant details.

## raw, but processed, data

### [gctree-node-data.csv](https://github.com/matsengrp/gcreplay/blob/main/results/gctree-node-data.csv)
Contains detailed information about each node in the phylogenetic trees inferred for each family tree. query the "naive_reversions_first" column values to get the data used for the primary analysis in the manuscript.

### [observed-seqs.csv](https://github.com/matsengrp/gcreplay/blob/main/results/observed-seqs.csv)
A comprehensive collection of all B-cell sequences observed in the study.

## Raw (raw) Data

### [data/NGS-gz](https://github.com/matsengrp/gcreplay/blob/main/data/NGS-gz)
Directory containing compressed Next Generation Sequencing (NGS) data files that serve as the raw input for the main pipeline. These files contain the sequencing reads from B-cell receptor repertoires before processing and analysis.

## (Notebook) Analysis Results

### [results/notebooks/cell-summaries/naive_reversions_first/cell_table.csv](https://github.com/matsengrp/gcreplay/blob/main/results/notebooks/cell-summaries/naive_reversions_first/cell_table.csv)
<!-- TODO WILL -->

### [results/notebooks/cell-summaries/naive_reversions_first/gc_summary.csv](https://github.com/matsengrp/gcreplay/blob/main/results/notebooks/cell-summaries/naive_reversions_first/gc_summary.csv)
<!-- TODO WILL -->

### [results/notebooks/mutations/naive_reversions_first/data.csv](https://github.com/matsengrp/gcreplay/blob/main/results/notebooks/mutations/naive_reversions_first/data.csv)
<!-- TODO WILL -->

### [results/notebooks/NDS-LB/naive_reversions_first/data.csv](https://github.com/matsengrp/gcreplay/blob/main/results/notebooks/NDS-LB/naive_reversions_first/data.csv)
<!-- TODO WILL -->

## Phylogenetic Trees

### [results/gctrees](https://github.com/matsengrp/gcreplay/blob/main/results/gctrees)

Directory containing the phylogenetic trees generated from the B-cell receptor sequences, for each germinal center according to it's "uid" (see [the gc metadata file](https://github.com/matsengrp/gcreplay/blob/main/gc_metadata.csv)).
In the manuscript primary analysis, we used the files under the /naive_reversions_first/ directories e.g. for tree `D15_M10_GC23` (Day 15, Mouse 10, Germinal Center 23), the files are:

```
results/gctrees/D15_M10_GC23
├── default                             <--- supplementary analysis
│   ├── gctree.inference.1.nk
│   ├── gctree.inference.1.svg
│   ├── gctree.p
│   └── mut_seq_annotated_nodes.svg
├── naive_reversions_first              < ---- * Primary analysis *
│   ├── gctree.inference.1.nk
│   ├── gctree.inference.1.svg
│   ├── gctree.p
│   └── mut_seq_annotated_nodes.svg
└── naive_reversions_no_bp              <--- supplementary analysis
    ├── gctree.inference.1.nk
    ├── gctree.inference.1.svg
    ├── gctree.p
    └── mut_seq_annotated_nodes.svg
```

 `gctree.p` contains the phylogenetic data in a pickle format that can be read by [gctree](https://github.com/matsengrp/gctree) version 4.3.0, The other files are a newick formatted tree, and two .svg visualizations show the annotated nodes of the tree including mutations, and a sequence id that can be found within the [gctree-node-data.csv](https://github.com/matsengrp/gcreplay/blob/main/results/gctree-node-data.csv) file.

## Accessing the Results

You could just use the GitHub web interface for accessing singular files of interest, but  here are further instructions for the detail oriented.

Many of the result files are stored using Git Large File Storage (LFS) due to their size. To access these files, you'll need to [install Git LFS](https://git-lfs.com/) first and then clone the repository:

```bash
git clone https://github.com/matsengrp/gcreplay.git
cd gcreplay
```

This will automatically pull the LFS files along with the regular repository contents. If you had already cloned the repository before installing Git LFS, you can pull the LFS files with:
```bash
git lfs pull
```