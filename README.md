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
    - In fact it might have increased B cells with affinity increases (more efficient affinity maturation?) but don’t remember if there was statistical significance