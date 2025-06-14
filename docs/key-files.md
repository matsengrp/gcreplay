# GCReplay Pipeline Results Documentation

This document provides an overview of the key files in the [gcreplay repository](https://github.com/matsengrp/gcreplay).

## Accessing the results

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

## Manuscript Figures:

Here, we tie the figures seen in our manuscript to the source code/data from which they were derived.

### Chapter 1
- Figures 1(C-J), S1(E), and S2(A-B) were generated from [NDS-LB.ipynb](https://github.com/matsengrp/gcreplay/blob/main/analysis/NDS-LB.ipynb).

### Chapter 2
- Figures 2(A,B,H,I), and S3(C) were generated from [mutations.ipynb](https://github.com/matsengrp/gcreplay/blob/main/analysis/mutations.ipynb).
- Figure S3(J) was generated from [NDS-LB.ipynb](https://github.com/matsengrp/gcreplay/blob/main/analysis/NDS-LB.ipynb).
- Figure S3(B) was generated from [collapse_scores.mb](https://github.com/jbloomlab/Ab-CGGnaive_DMS/blob/main/results/summary/collapse_scores.md).

### Chapter 3
- Figure 3(A) was generated from the analysis in [passenger](https://github.com/matsengrp/gcreplay/tree/main/passenger), with final plots in [plots_for_paper.ipynb](https://github.com/matsengrp/gcreplay/blob/main/passenger/plots_for_paper.ipynb).
- Figures 3(B-I), and S4(C-D) were generated from analysis in [mutations.ipynb](https://github.com/matsengrp/gcreplay/blob/main/analysis/mutations.ipynb).
- Figures 3(B,F) was generated from interactive software using the key files [here](https://matsen.group/gcreplay/key-files/#raw-processed-data).
- Figure 3(G) was generated from analysis in [mutation-accessibility.ipynb](https://github.com/matsengrp/gcreplay/blob/main/analysis/mutation-accessibility.ipynb).
- Figures 3(H-I) were generated from interactive software using the key files [here](https://matsen.group/gcreplay/key-files/#raw-processed-data).

### Chapter 4
- Figure 4(A) was generated from the analysis in [cell-summaries.ipynb](https://github.com/matsengrp/gcreplay/blob/main/analysis/cell-summaries.ipynb).
- Figures 4(B-E) were generated from analysis in [fitness-trajectories](https://github.com/matsengrp/gcreplay/blob/main/analysis/phenotype-trajectories.ipynb).
- Figures 5(A,C,F,G) were generated from analysis in [NDS-LB.ipynb](https://github.com/matsengrp/gcreplay/blob/main/analysis/NDS-LB.ipynb).
- Figures 5(B,D,H) were generated from analysis in [fitness-regression.ipynb](https://github.com/matsengrp/gcreplay/blob/main/analysis/fitness-regression.ipynb).
- Figure 5(E) was generated from interactive software using the key files [here](https://matsen.group/gcreplay/key-files/#raw-processed-data).
- Figure S5(A) was generated from [cell-summaries.ipynb](https://github.com/matsengrp/gcreplay/blob/main/analysis/cell-summaries.ipynb).
- Figure S5(C) was generated from interactive software using the key files [here](https://matsen.group/gcreplay/key-files/#raw-processed-data).
- Figure S5(B) was generated from [NDS-LB.ipynb](https://github.com/matsengrp/gcreplay/blob/main/analysis/NDS-LB.ipynb).
- Figure S5(D) was generated from [fitness-regression.ipynb](https://github.com/matsengrp/gcreplay/blob/main/analysis/fitness-regression.ipynb).

### Chapter 5
- Figure 6(A-H) were all generated from [affinity-fitness-response.ipynb](https://github.com/matsengrp/gcreplay/blob/main/analysis/affinity-fitness-response.ipynb)


## Sample Metadata 

### [gc_metadata.csv](https://github.com/matsengrp/gcreplay/blob/main/gc_metadata.csv)

The metadata for germinal centers analyzed across all replay experiments is stored in this file. The columns are as follows:

- **uid**: (_unique identifier_) This is in the format of D<imm_duration>_M<mouse>_GC<gc>. 
- **ngs_id**: (_ngs file id_) reference to associated sequencing run specified in [ngs_manifest.csv](https://github.com/matsengrp/gcreplay/blob/main/ngs_manifest.csv).
- **imm_duration**: (_immunization duration_) This is the number of weeks after immunization that the sample was taken.
- **mouse**: (_mouse number_) This is the unique identifier for the mouse.
- **gc**: (_germinal center number_) This is the unique identifier for the germinal center.
- **strain**: (_mouse strain_) This is the strain of the mouse.
- **node**: (_lymph node number_) This identifies the specific lymph node location in the mouse.
- **cell_type**: (_cell type_) This is the type of cell that was sequenced.
- **plate**: (_plate barcode_) This is the barcode of the plate that the sample was sequenced on.
- **hc_barcode** : (_heavy chain barcode_) This is the barcode of the heavy chain.
- **lc_barcode** : (_light chain barcode_) This is the barcode of the light chain.
- **row** : (_well rows_) Which rows of the 96-well plate contain the sample.
- **col** : (_well columns_) Which columns of the 96-well plate contain the sample.


### [ngs_manifest.csv](https://github.com/matsengrp/gcreplay/blob/main/ngs_manifest.csv)
Serves as a manifest for the Next Generation Sequencing (NGS) data, detailing the sequencing runs and their corresponding metadata.

## Raw (processed) data

### [gctree node data table](https://github.com/matsengrp/gcreplay/blob/main/results/gctree-node-data.csv)
Contains detailed information about each node in the phylogenetic trees inferred for each family tree. query the "naive_reversions_first" column values to get the data used for the primary analysis in the manuscript.

### [Observed BCR's table](https://github.com/matsengrp/gcreplay/blob/main/results/observed-seqs.csv)
A comprehensive collection of all B-cell sequences observed in the study.
This includes a ton of information about each sequence, including the sequence itself, the number of times it was observed, the mutations it has aquired, and the germinal center it was observed in.

### [Phylogenetic `gctree` directories](https://github.com/matsengrp/gcreplay/blob/main/results/gctrees)

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


## (Notebook) Analysis Results

### [results/notebooks/cell-summaries/naive_reversions_first/cell_table.csv](https://github.com/matsengrp/gcreplay/blob/main/results/notebooks/cell-summaries/naive_reversions_first/cell_table.csv)
A comprehensive dataset containing single-cell information from all germinal centers (GCs) analyzed in the study. Each row represents a single B cell with detailed sequence data, including unique identifiers, mutation information (nucleotide and amino acid substitutions), antigen binding affinity (delta_bind), gene expression measurements (delta_expr), and immunoglobulin isotype information. This file serves as the foundation for cell-level analyses throughout the study.

### [results/notebooks/cell-summaries/naive_reversions_first/gc_summary.csv](https://github.com/matsengrp/gcreplay/blob/main/results/notebooks/cell-summaries/naive_reversions_first/gc_summary.csv)
An aggregated summary of germinal center characteristics, with each row representing a distinct germinal center. Contains metrics such as sampling time, number of cells sampled, median affinity measurements, dominance scores (indicating clonal diversity), and maximum Revving Evolutionary Influence (REI) scores. This file provides a high-level view of evolutionary dynamics and selection pressures within each germinal center structure.

### [results/notebooks/mutations/naive_reversions_first/data.csv](https://github.com/matsengrp/gcreplay/blob/main/results/notebooks/mutations/naive_reversions_first/data.csv)
A detailed catalog of all mutations observed across B-cell receptors, including their effects on antigen binding and gene expression. Each row represents a unique mutation with information about its location (position, chain), frequency (events, abundance), biophysical properties (binding affinity, expression changes), and spatial characteristics (distance to antigen). This dataset is crucial for understanding mutation patterns and selection pressures during B-cell affinity maturation.

### [results/notebooks/NDS-LB/naive_reversions_first/data.csv](https://github.com/matsengrp/gcreplay/blob/main/results/notebooks/NDS-LB/naive_reversions_first/data.csv)
Contains tree shape statistics and evolutionary metrics for each germinal center, generated from phylogenetic analysis. The data includes Normalized Dominance Scores (NDS) that quantify clonal diversity, Revving Evolutionary Influence (REI) measurements that identify key evolutionary events, and various affinity metrics that track selection for improved antigen binding. This file is essential for understanding the relationship between tree topology and functional outcomes during B-cell evolution.



## Raw (NGS) Data

### [NGS-gz](https://github.com/matsengrp/gcreplay/blob/main/data/NGS-gz)
Directory containing compressed Next Generation Sequencing (NGS) data files that serve as the raw input for the main pipeline. These files contain the sequencing reads from B-cell receptor repertoires before processing and analysis.
