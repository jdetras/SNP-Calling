#!/usr/bin/env nextflow
/*
 * SNP Calling Pipeline
 * Pipeline based on GATK best practices
 *
 * Jeffrey Detras
 */

nextflow.enable.dsl = 2

/*
 * Define default parameters
 */

params.genome = "$projectDir/data/reference/genome.fa"
params.reference = "$projectDir/data/reference/"
params.reads = "$projectDir/data/reads/palay/*-?_{1,2}.fq.gz"
params.output = "$projectDir/output"

genome_fai = file(params.genome + ".fai")
genome_bwt = file(params.genome + ".bwt")
genome_amb = file(params.genome + ".amb")
genome_ann = file(params.genome + ".ann")
genome_pac = file(params.genome + ".pac")
genome_sa  = file(params.genome + ".sa")
genome_dic = file(params.reference + "/genome.dict")

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
    alignment;
    sort_sam;
    fix_mate;
    mark_duplicates;
    add_read_groups;
    merge_bam;
    haplotype_caller;
} from './modules.nf'

/*
 * main pipeline
 */

workflow{
        
    //read alignment to reference using BWA 
    reads_ch = Channel.fromFilePairs(params.reads)
    
    alignment(
        params.genome,
        reads_ch,
        genome_amb, 
        genome_ann, 
        genome_bwt, 
        genome_pac, 
        genome_sa 
    )

    //processing of BAM files using PicardTools
    
    sort_sam(alignment.out)
    fix_mate(sort_sam.out)
    mark_duplicates(fix_mate.out)
    add_read_groups(mark_duplicates.out)

    //Merging multiple sample files using Samtools
    
    merge_bam(params.output,add_read_groups.out)

    //Variant calling using GATK

    haplotype_caller(
        params.genome, 
        genome_fai, 
        genome_dic, 
        merge_bam.out
    )
}
