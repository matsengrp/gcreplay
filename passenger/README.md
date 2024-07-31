# Passenger allele analysis

## Installation

Here's a working set of steps for a conda environment:

```
mamba create -n gcreplay python=3.9
conda activate gcreplay
pip install seaborn biopython pandas ipykernel scipy
conda install -c bioconda blast
```

## To replicate analysis

All instructions are relative to `gcreplay/passenger/`.

Download these files into `input`:

    outname_5_IgG_S4_R1_001_atleast-2.fastq.gz
    outname_5_IgM_S1_R1_001_atleast-2.fastq.gz
    outname_6_IgG_S5_R1_001_atleast-2.fastq.gz
    outname_6_IgM_S2_R1_001_atleast-2.fastq.gz
    outname_7_IgG_S6_R1_001_atleast-2.fastq.gz
    outname_7_IgM_S3_R1_001_atleast-2.fastq.gz
    outname_5_S1_R1_001_atleast-2.fastq.gz
    outname_6A_S2_R1_001_atleast-2.fastq.gz
    outname_7B_S3_R1_001_atleast-2.fastq.gz
    outname_AV_IgYKSTOP_11A_S1_R1_001_atleast-2.fastq.gz
    outname_AV_IgYKSTOP_8A_S2_R1_001_atleast-2.fastq.gz
    outname_AV_IgYKSTOP_9A_S3_R1_001_atleast-2.fastq.gz


Enter the `passenger-blast` directory and execute

    bash process-and-blast.sh

Execute the cells in these notebooks to do processing, as well as some plotting.

* `passenger-8a.ipynb`
* `passenger-9a.ipynb`
* `passenger-11a.ipynb`
* `passenger-igh.ipynb`

Then execute these notebooks, which perform model inference across sequencing runs.

* `igh_passenger_aggregate.ipynb`
* `igk_passenger_aggregate.ipynb`

The final plots for the paper are made in `plots_for_paper.ipynb`.


## Key directories
* `input`: where input files should go to replicate analysis
* `output`: output files
* `passenger-blast`: where we use BLAST to identify sequences that come from passenger allele