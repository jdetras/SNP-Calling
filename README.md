# SNP-Calling
GATK Variant calling pipeline for genomic data using Nextflow

[![nextflow](https://img.shields.io/badge/nextflow-%E2%89%A522.04.5-brightgreen.svg)](http://nextflow.io)

## Quickstart

Install Nextflow using the following command: 

  `curl -s https://get.nextflow.io | bash`
  
Index reference genome

  `$ bwa index /path/to/reference/genome.fa`
 
  `$ samtools faidx /path/to/reference/genome.fa`
  
  `$ gatk CreateSequenceDictionary -R /path/to/genome.fa -O genome.dict`

Launch the pipeline execution with the following command:

  `nextflow run jdetras/snp-calling -profile docker`
  
## Pipeline Description

The variant calling pipeline follows the recommended practices from GATK. The input genomic data are aligned to a reference genome using BWA. The alignemnt files are processed using Picard Tools. Variant calling is done using samtools and GATK. 

## Input files

The input files required to run the pipeline:
* Genomic sequence paired reads, `*_{1,2}.fq.gz`
* Reference genome, `*.fa`

## Pipeline parameters

### Usage
Usage: `nextflow run jdetras/snp-calling -profile docker [options]`

Options:

* `--reads` 
* `--genome`
* `--output`

Example: 
  `$ nextflow run jdetras/snp-calling -profile docker --reads '/path/to/reads/*_{1,2}.fq.gz' --genome '/path/to/reference/genome.fa' --output '/path/to/output'`

#### `--reads`

* The path to the FASTQ read files.
* Wildcards (*, ?) can be used to declare multiple reads. Use single quotes when wildcards are used. 
* Default parameter: `$projectDir/data/reads/*_{1,2}.fq.gz`

Example: 
  `$ nextflow run jdetras/snp-calling -profile docker --reads '/path/to/reads/*_{1,2}.fq.gz'`
  
#### `--genome`

* The path to the genome file in fasta format.
* The extension is `.fa`.
* Default parameter: `$projectDir/data/reference/genome.fa`

Example:
  `$ nextflow run jdetras/snp-calling -profile docker --genome /path/to/reference/genome.fa`
    
#### `--output`

* The path to the directory for the output files.
* Default parameter: `$projectDir/output`

## Software

* [BWA 0.7.17](http://bio-bwa.sourceforge.net/)
* [Samtools 1.3.1](http://www.htslib.org/)
* [GATK 4.2.6.1](https://gatk.broadinstitute.org/) 