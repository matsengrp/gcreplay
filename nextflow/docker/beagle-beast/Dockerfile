
FROM quay.io/hdc-workflows/ubuntu:20.04

LABEL maintainer "Jared Galloway <jgallowa@fredhutch.rg>" \
      version "1.0" \
      description "GC Replay Analysis Workflow"

# install needed tools
RUN apt-get update --fix-missing -qq && \
    DEBIAN_FRONTEND=noninteractive \
    apt-get install -y -q \
    cmake \
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

# install python deps
RUN apt-get install -y wget && rm -rf /var/lib/apt/lists/*
RUN wget \
    https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    && mkdir /opt/.conda \
    && bash Miniconda3-latest-Linux-x86_64.sh -b -p /opt/miniconda3 \
    && rm -f Miniconda3-latest-Linux-x86_64.sh
    # ./Miniconda3-py38_4.8.3-Linux-x86_64.sh -b\

# install phylip
RUN conda install pip
RUN pip install biopython==1.81

RUN git clone --depth=1 https://github.com/beagle-dev/beagle-lib.git \
    && mkdir beagle-lib/build \
    && (cd beagle-lib/build && cmake -DBUILD_OPENCL=OFF -DCMAKE_INSTALL_PREFIX:PATH=/ .. && make install)

ENV LD_LIBRARY_PATH=/lib
ENV PKG_CONFIG_PATH=/lib/pkgconfig

# install BEAST
RUN curl -fsSL --proto '=https' \
    https://github.com/beast-dev/beast-mcmc/releases/download/v1.10.4/BEASTv1.10.4.tgz \
    --output BEASTv1.10.4.tgz \
    && tar -zxvf BEASTv1.10.4.tgz \
    && (cd /bin && ln -s /BEASTv1.10.4/bin/beast ./beast)
    
# install BEASTGEN
# Note: this was downloaded from google drive which doesn't easily allow 
# for curling a specific file. The command below returns a 400 error
#RUN curl -fsSL --proto '=https' 
#    https://drive.google.com/file/d/0B37cqWL7UhTAVFVhQ2o1Y093b1k \
#    --output BEASTGen_v1.0.2.tgz \ 
#    && tar -xvf BEASTGen_v1.0.2.tgz \
#    && (cd /bin && ln -s /BEASTGen_v1.0.2/bin/beastgen ./beastgen)

COPY ./BEASTGen_v1.0.2.tgz .
RUN tar -zxvf BEASTGen_v1.0.2.tgz \
    && chmod +x /BEASTGen_v1.0.2/bin/beastgen \
    && (cd /bin && ln -s /BEASTGen_v1.0.2/bin/beastgen ./beastgen)
