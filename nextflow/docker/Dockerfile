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
    pkg-config \
    python3 \
    python3-pip \ 
    python3-venv

# create a virtual environment for python tools
ENV VIRTUAL_ENV=/opt/venv
RUN python3 -m venv $VIRTUAL_ENV
ENV PATH="$VIRTUAL_ENV/bin:$PATH"

# install cutadapt
RUN python3 -m pip install --upgrade cutadapt

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



