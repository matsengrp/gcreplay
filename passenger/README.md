# Passenger allele analysis

## To replicate analysis

All instructions are relative to `gcreplay/passenger/`.

Download these files into `input`:

    5_IgG_S4_R1_001_atleast-2.fastq.gz
    5_IgM_S1_R1_001_atleast-2.fastq.gz
    6_IgG_S5_R1_001_atleast-2.fastq.gz
    6_IgM_S2_R1_001_atleast-2.fastq.gz
    7_IgG_S6_R1_001_atleast-2.fastq.gz
    7_IgM_S3_R1_001_atleast-2.fastq.gz
    AV_IgYKSTOP_11A_collapse-unique.fastq.gz
    AV_IgYKSTOP_8A_collapse-unique.fastq.gz
    AV_IgYKSTOP_9A_collapse-unique.fastq.gz

Enter the `passenger-blast` directory and execute

    bash process-and-blast.sh

Execute the cells in these notebooks to do processing, as well as some plotting.

* `passenger-8a.ipynb`
* `passenger-9a.ipynb`
* `passenger-11a.ipynb`
* `passenger-igh.ipynb`

Then execute these notebooks to make final plots.

* `igh_passenger_aggregate.ipynb`
* `igk_passenger_aggregate.ipynb`


## Key directories
* `input`: where input files should go to replicate analysis
* `output`: output files
* `passenger-blast`: where we use BLAST to identify sequences that come from passenger allele