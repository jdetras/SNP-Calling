#!/usr/bin/env nextflow
/*
 * SNP Calling Pipeline
 * 
 * Pipeline based on GATK best practices
 *
 * Jeffrey Detras
 */

nextflow.enable.dsl = 2
/*
 * Define default parameters
 */
params.genome = "$baseDir/data/reference/Chr7.fasta"
params.reference = "$baseDir/data/reference"
params.reads = "$baseDir/data/reads/palay/*-?_{1,2}.fq.gz"
params.data = "$baseDir/data"
params.output = "$baseDir/data/output"

genome_fai = file(params.genome + ".fai")
genome_bwt = file(params.genome + ".bwt")
genome_amb = file(params.genome + ".amb")
genome_ann = file(params.genome + ".ann")
genome_pac = file(params.genome + ".pac")
genome_sa  = file(params.genome + ".sa")

log.info """\
SNP CALLING
===========
genome : $params.genome
reads  : $params.reads
output : $params.output
"""

/*
 * Import modules
 */
include {
//	prepare_reference;
	alignment;
	sort_sam;
	fix_mate;
	mark_duplicates;
	add_read_groups;
//	merge_bam;
//	haplotype_caller;
} from './modules.nf'

/*
 * main
 */
workflow{
	reads_ch = Channel.fromFilePairs(params.reads)
	bams = Channel.fromPath("$params.output/*.rg.bam")
//	prepare_reference(params.genome)
        alignment(
		params.genome,
		reads_ch,
		genome_amb, 
		genome_ann, 
		genome_bwt, 
		genome_pac, 
		genome_sa )
	sort_sam(alignment.out )
	fix_mate(sort_sam.out )
	mark_duplicates(fix_mate.out)
	add_read_groups(mark_duplicates.out)
//	merge_bam(bams)
//	haplotype_caller
}
