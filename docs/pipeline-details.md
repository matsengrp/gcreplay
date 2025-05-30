
# Replay analysis pipeline

We provide a Nextflow pipeline for running the vast majority of the analyses seen in our manuscript. 
There are a few analyses that were run outside of this pipeline including
the [passenger analysis](https://github.com/matsengrp/gcreplay/tree/main/passenger), and
the [POTP analysis](https://github.com/matsengrp/gcreplay/blob/main/analysis/affinity-fitness-response.ipynb).


## Reproduce analysis

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
    $ nextflow run main.nf -profile docker -resume

Note that this pipeline is computationally intensive, we run the pipeline on a SLURM cluster using the configurations seen in 
[nextflow.config](https://github.com/matsengrp/gcreplay/blob/main/nextflow.config)


## Pipeline parameters

* `--ngs_manifest` - Path to the NGS manifest CSV file containing sequencing run information. Default: `ngs_manifest.csv`

* `--gc_metadata` - Path to the germinal center metadata CSV file. Default: `gc_metadata.csv`

* `--reads_prefix` - Directory containing the NGS read files. Default: `data/NGS-gz`

* `--results` - Directory where pipeline results will be stored. Default: `results/`

* `--plate_barcodes` - File containing plate barcode sequences. Default: `data/barcodes/plateBC.txt`

* `--well_barcodes` - File containing 96-well plate barcode sequences. Default: `data/barcodes/96FBC.txt`

* `--partis_anno_dir` - Directory containing partis annotation germline data. Default: `data/partis_annotation/germlines`

* `--hdag_sub` - File containing HDAG substitution model data. Default: `data/mutability/MK_RS5NF_substitution.csv`

* `--hdag_mut` - File containing HDAG mutability model data. Default: `data/mutability/MK_RS5NF_mutability.csv`

* `--chigy_hc_mut_rates` - File containing chimeric gamma heavy chain mutation rates. Default: `data/mutability/chigy_hc_mutation_rates_nt.csv`

* `--chigy_lc_mut_rates` - File containing chimeric gamma light chain mutation rates. Default: `data/mutability/chigy_lc_mutation_rates_nt.csv`

* `--pdb` - PDB structure file for antibody analysis. Default: `data/AbCGG_structure/combined_ch2_eh2-coot_IMGT.pdb`

* `--dms_vscores` - File containing deep mutational scanning variant scores. Default: `data/dms/final_variant_scores.csv`

* `--dms_sites` - File containing naive sites information. Default: `data/dms/CGGnaive_sites.csv`

* `--heavy_chain_motif` - DNA motif sequence used to identify heavy chain reads. Default: `aGCgACgGGaGTtCAcagACTGCAACCGGTGTACATTCC`

* `--light_chain_motif` - DNA motif sequence used to identify light chain reads. Default: `aGCgACgGGaGTtCAcagGTATACATGTTGCTGTGGTTGTCTG`

* `--igk_idx` - Index position for immunoglobulin kappa chain processing. Default: `336`

* `--bcr_count_thresh` - Minimum count threshold for B-cell receptor sequences. Default: `5`

## Pipeline steps

(1) trim and combine paired end files
(2) demultiplex both plates _and_ wells
(3) Split Heavy and light chain reads per well
(4) Collapse identical sequences, while retaining and ordering rank of each sequence per well
(5) Prune to keep the top N sequences observed from each well
(6) Merge the results, formatting for [partis](https://github.com/psathyrella/partis) annotation.
(7) [partis](https://github.com/psathyrella/partis) annotation
(8) curation, cleaning, and merging of Heavy and Light Chains
(9) [gctree](https://github.com/matsengrp/gctree) lineage inference using HDAG.
(10) Merge the results

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