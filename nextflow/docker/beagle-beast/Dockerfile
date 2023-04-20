
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

RUN git clone --depth=1 https://github.com/beagle-dev/beagle-lib.git \
    && mkdir beagle-lib/build \
    && (cd beagle-lib/build && cmake -DBUILD_OPENCL=OFF -DCMAKE_INSTALL_PREFIX:PATH=$HOME .. && make install)

ENV LD_LIBRARY_PATH=$HOME/lib:$LD_LIBRARY_PATH
ENV PKG_CONFIG_PATH=$HOME/lib/pkgconfig:$PKG_CONFIG_PATH

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
    && (cd /bin && ln -s /BEASTGen_v1.0.2/bin/beastgen ./beastgen)
