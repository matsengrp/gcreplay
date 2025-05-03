# GC replay

This repository contains the code and data for the analysis of B cell receptor (BCR) sequences isolated in Germinal Center (GC) - Replay experiments.

For more details, please see the analysis [website](https://matsengrp.github.io/gcreplay/).

## Reproduce analysis

Install `Nextflow` by using the following command:

    $ curl -s https://get.nextflow.io | bash

Download the `Docker` Desktop, there exists several distributions packaged for
various linux flavors

    $ curl -fsSL https://get.docker.com -o get-docker.sh && sudo sh get-docker.sh

Note: the [Dockerfile](docker/Dockerfile) contains all the required dependencies.
Add the `-profile docker` to enable the containerized execution to the
example command line shown below.

Launch the pipeline execution with the following command:

    $ git clone git@github.com:matsengrp/gcreplay.git && cd gcreplay
    $ nextflow run main.nf -profile docker -resume

## GC Metadata

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

## Building docs

We use `pixi` to manage dependecies. To build the documentation, first install `pixi`:

```bash
$ curl -fsSL https://pixi.sh/install.sh | bash
$ cd gcreplay
$ pixi install
```

Then, run the following command to build the documentation:

```bash
mkdocs serve
```