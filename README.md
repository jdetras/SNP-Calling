# SNP-Calling
GATK Variant calling pipeline for genomic data using Nextflow

[![nextflow](https://img.shields.io/badge/nextflow-%E2%89%A522.04.5.5708-brightgreen.svg)](http://nextflow.io)

## Install Nextflow using the following command: 

  `curl -s https://get.nextflow.io | bash`
  
## Pipeline Description

The variant calling pipeline follows the recommended practices from GATK. The input genomic data are aligned to a reference genome using BWA. The alignemnt files are processed using Picard Tools. Variant calling is done using samtools and GATK. 

## Input files

The input files required to run the pipeline:
* Genomic sequence paired reads, `*_{1,2}.fq.gz`
* Reference genome, `*.fa`

## Pipeline parameters

### Usage
Usage: `nextflow run main.nf [options]`

Options:

* `--reads` 
* `--genome`
* `--output`

Example: 
`$ nextflow run main.nf --reads '/home/data/reads/*_{1,2}.fq.gz' --genome '/home/data/reference/genome.fa' --output '/home/output'`

Note: `main.nf` requires `module.nf` on the same working directory

#### `--reads`

* The path to the FASTQ read files.
* Wildcards (*, ?) can be used to declare multiple reads. Use single quotes when wildcards are used. 
* Default parameter: `$projectDir/data/reads/*_{1,2}.fq.gz`

Example: 
  `$ nextflow run main.nf --reads '/home/data/reads/*_{1,2}.fq.gz'`
  
#### `--genome`

* The path to the genome file in fasta format.
* The extension is `.fa`.
* Default parameter: `$projectDir/data/reference/genome.fa`

Example:
  `$ nextflow run main.nf --genome home/data/reference/genome.fa`
  
Note: indexing of reference genome should be done before running the Nextflow pipeline

  `$ bwa index /home/data/reference/genome.fa`
 
  `$ samtools faidx /home/data/reference/genome.fa`
  
  `$ picard CreateSequenceDictionary -R genome.fa -O genome.dict`
  
#### `--output`

* The path to the directory for the output files.
* Default parameter: `$projectDir/output`

## Software

* [BWA](http://bio-bwa.sourceforge.net/)
* [Picard Tools](https://broadinstitute.github.io/picard/)
* [GATK](https://gatk.broadinstitute.org/)
* [Samtools](http://www.htslib.org/)
