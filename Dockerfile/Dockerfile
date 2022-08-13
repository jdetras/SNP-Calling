# get java 
FROM openjdk:latest as JavaBuilder

# second stage
FROM debian:stable

LABEL maintainer "Emilio Palumbo <emilio.palumbo@crg.eu>" \
      version "1.0" \
      description "Varian Calling Analysis with RNA-seq data"

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

# install VCFtools
RUN curl -fksSL https://github.com/vcftools/vcftools/releases/download/v0.1.14/vcftools-0.1.14.tar.gz | tar xz && \
    cd vcftools-0.1.14 && \
    ./configure && make && make install

RUN apt-get update && apt-get install -y python \
  && curl -fsSL https://github.com/broadinstitute/gatk/releases/download/4.2.6.1/gatk-4.2.6.1.zip > gatk-4.zip \
  && unzip gatk-4.zip \
  && rm gatk-4.zip

RUN apt-get install bwa

COPY --from=JavaBuilder /usr/java/openjdk-18 /usr/java/openjdk-18

ENV JAVA_HOME=/usr/java/openjdk-18

ENV PATH=/gatk-4.2.6.1:/usr/java/openjdk-18/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin

ENV LANG=C.UTF-8

ENV JAVA_VERSION=18.0.2