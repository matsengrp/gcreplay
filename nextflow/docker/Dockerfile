FROM quay.io/hdc-workflows/ubuntu:20.04

LABEL maintainer "Jared Galloway <jgallowa@fredhutch.rg>" \
      version "1.0" \
      description "GC Replay Analysis Workflow"

# install needed tools
RUN apt-get update --fix-missing -qq && \
    DEBIAN_FRONTEND=noninteractive \
    apt-get install -y -q \
    git \
    curl \
    locales \
    libncurses5-dev  \
    libncursesw5-dev \
    build-essential \
    libtool \
    automake \
    zlib1g-dev \
    libbz2-dev \
    pkg-config

# install graphical libraries used by qt and vispy
RUN apt-get install -y xvfb
RUN apt-get install -qqy mesa-utils libgl1-mesa-glx  libglib2.0-0
RUN apt-get install -qqy libfontconfig1 libxrender1 libdbus-1-3 libxkbcommon-x11-0 libxi6
RUN apt-get install -qqy libxcb-icccm4 libxcb-image0 libxcb-keysyms1 libxcb-randr0 libxcb-render-util0
RUN apt-get install -qqy libxcb-xinerama0 libxcb-xinput0 libxcb-xfixes0 libxcb-shape0

# For this to work with a singularity wrapper
# we have top move it somewhere other than root
ENV PATH="/opt/miniconda3/bin:${PATH}"
ARG PATH="/opt/miniconda3/bin:${PATH}"

RUN apt-get install -y wget && rm -rf /var/lib/apt/lists/*
RUN wget \
    https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    && mkdir /opt/.conda \
    && bash Miniconda3-latest-Linux-x86_64.sh -b \
    ./Miniconda3-py38_4.8.3-Linux-x86_64.sh -b -p /opt/miniconda3 \
    && rm -f Miniconda3-latest-Linux-x86_64.sh

# install pip
RUN conda install pip

# install phylip
RUN conda install -c bioconda phylip

# Other pip deps
RUN pip install --upgrade \
    cutadapt \
    regex \
    click \
    biopython

# install branch specific version of gctree
RUN pip install git+https://github.com/matsengrp/gctree.git@fit-on-history-dag

# install libgtextutils
RUN curl -L 'https://github.com/agordon/libgtextutils/releases/download/0.7/libgtextutils-0.7.tar.gz' \
  --output libgtextutils-0.7.tar.gz \
  && gzip -d libgtextutils-0.7.tar.gz \
  && tar -xvf libgtextutils-0.7.tar \
  && (cd libgtextutils-0.7 && ./configure && make CXXFLAGS='-std=c++03 -O1') \
  && (cd libgtextutils-0.7 && make install)

# install fastx toolkit
RUN curl -L 'https://github.com/agordon/fastx_toolkit/releases/download/0.0.14/fastx_toolkit-0.0.14.tar.bz2' \
  --output fastx_toolkit-0.0.14.tar.bz2 \
  && tar -xjvf fastx_toolkit-0.0.14.tar.bz2 \
  && (cd fastx_toolkit-0.0.14 && ./configure && make CXXFLAGS="-Wno-error") \
  && (cd fastx_toolkit-0.0.14 && make install)

# install pandaseq
RUN git clone https://github.com/neufeld/pandaseq.git \
  && (cd pandaseq && ./autogen.sh && ./configure && make && make install && ldconfig)
