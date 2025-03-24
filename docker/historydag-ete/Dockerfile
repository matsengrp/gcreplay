
FROM quay.io/hdc-workflows/ubuntu:20.04

LABEL maintainer "Jared Galloway <jgallowa@fredhutch.rg>" \
      version "1.0" \
      description "GC Replay Analysis Workflow"

# install needed tools
RUN apt-get update --fix-missing -qq && \
    DEBIAN_FRONTEND=noninteractive \
    apt-get install -y -q \
    build-essential \
    autoconf \
    automake \
    libtool \
    git \
    curl \
    pkg-config \
    openjdk-11-jdk \
    openjdk-8-jre \
    openjdk-8-jdk

# For this to work with a singularity wrapper
# we have top move it somewhere other than root
ENV PATH="/opt/miniconda3/bin:${PATH}"
ARG PATH="/opt/miniconda3/bin:${PATH}"

# wget for anaconda
RUN apt-get install -y wget && rm -rf /var/lib/apt/lists/*

# install to the directory set above
RUN wget \
    https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    && mkdir /opt/.conda \
    && bash Miniconda3-latest-Linux-x86_64.sh -b -p /opt/miniconda3 \
    && rm -f Miniconda3-latest-Linux-x86_64.sh

# TODO if you want viz, we need to get python3.9 (<3.10)?
RUN conda install pip
RUN pip install pandas numpy biopython seaborn matplotlib ete3 dendropy regex

# currently the ete code lives only in a branch
RUN git clone https://github.com/matsengrp/historydag.git \
    && (cd historydag && git checkout wd-load-beast-trees && pip install -e .)
