# GC replay

This repository contains the code and data for the analysis of B cell receptor (BCR) sequences from germinal center (GC) B cells in mice.

For more details, please see the analysis [website](https://matsengrp.github.io/gcreplay/).

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
    $ nextflow run nextflow/main.nf -profile docker -resume