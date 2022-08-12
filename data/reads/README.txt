README

These are FASTQ files derived from sequencing a single sample ( IRIS-313_11573 )
For training purposes, these represent only a small fraction of the original reads.


The reads are separated into two "read groups". Sequencing reads coming e.g from different libraries should be analysed separately, and will be assigned different "read groups".
Each read group has its own pair of FASTQ files (r1,r2) and, according to best practices, should be processed separately (alignment and deduplication steps).  Before calling SNPs, the BAM files from different read groups will be merged into a single BAM file.

* * *

To know more about read groups, please refer to this document
https://gatk.broadinstitute.org/hc/en-us/articles/360035890671-Read-groups 

There is no formal definition of what a 'read group' is, however in practice this term refers to a set of reads that are generated from a single run of a sequencing instrument.

In the simple case where a single library preparation derived from a single biological sample was run on a single lane of a flow cell, all the reads from that lane run belong to the same read group. When multiplexing is involved, then each subset of reads originating from a separate library run on that lane will constitute a separate read group
