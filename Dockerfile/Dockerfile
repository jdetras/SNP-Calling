FROM debian:stable

LABEL maintainer "Jeffrey Detras" \
      version "1.0" \
      description "SNP Calling on genomic data using Nextflow"

# install needed tools
RUN apt-get update --fix-missing -qq && apt-get install -y -q \
    curl \
    locales \
    libncurses5-dev  \
    libncursesw5-dev \
    build-essential \
    pkg-config \
    zlib1g-dev \
    bzip2 \
    r-base \
    && apt-get clean \
    && apt-get purge \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*


# install SAMtools
RUN curl -fksSL https://github.com/samtools/samtools/releases/download/1.3.1/samtools-1.3.1.tar.bz2 | tar xj && \
    cd samtools-1.3.1 && \
    make all all-htslib && make install install-htslib

RUN apt-get update && apt-get install -y bwa

