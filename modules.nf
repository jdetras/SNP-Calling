process prepare_reference {
    tag "$genome.baseName"
    publishDir params.reference

    input:
      path genome
    
    output:
      tuple \
        path("${genome}.fai"), \
        path("${genome}.amb"), \
        path("${genome}.ann"), \
        path("${genome}.bwt"), \
        path("${genome}.pac"), \
        path("${genome}.sa")

    script:
    """
    samtools faidx $genome
    bwa index $genome
    """	
}

process alignment {
    publishDir params.output

    input:
      path genome
      tuple val(sampleName), file(reads)
      path genome_amb
      path genome_ann
      path genome_bwt
      path genome_pac
      path genome_sa

    output:
      file "${sampleName}.sam"

    script:
    """
    bwa mem -M ${genome} ${reads[0]} ${reads[1]} > ${sampleName}.sam
    """
}

process sort_sam {
    tag "$sam.baseName"
    publishDir params.output 

    input:
      path sam

    output:
      path "${sam.baseName}.sorted.bam"

    script:
    """
    picard SortSam \
      -I ${sam} \
      -O ${sam.baseName}.sorted.bam \
      -SO coordinate \
      --VALIDATION_STRINGENCY LENIENT \
      --CREATE_INDEX TRUE
    """
}

process fix_mate {
    tag "$bam.simpleName"
    publishDir params.output 
    
    input:
      path bam

    output:
      path "${bam.simpleName}.fixmate.bam"

    script:
    """
    picard FixMateInformation \
      -I ${bam} \
      -O ${bam.simpleName}.fixmate.bam \
      --VALIDATION_STRINGENCY LENIENT \
      --CREATE_INDEX TRUE
    """
}

process mark_duplicates {
    tag "$bam.simpleName"
    publishDir params.output 

    input: 
      path bam

    output:
      path "${bam.simpleName}.markdup.bam"

    script:
    """
    picard MarkDuplicates \
       -I ${bam} \
       -O ${bam.simpleName}.markdup.bam \
       --METRICS_FILE ${bam.simpleName}.metrics \
       --VALIDATION_STRINGENCY LENIENT \
       --CREATE_INDEX TRUE
    """
}

process add_read_groups {
    tag "$bam.simpleName"
    publishDir params.output

    input:
      path bam

    output:
      path "${bam.simpleName}.rg.bam"

    script:
    """
    picard AddOrReplaceReadGroups \
      -I ${bam} \
      -O ${bam.simpleName}.rg.bam \
      -ID ${bam.simpleName} \
      -PL Illumina \
      -SM rice \
      -LB 3k \
      -CN BGI \
      -PU gsl \
      --VALIDATION_STRINGENCY LENIENT \
      --CREATE_INDEX TRUE
    """
}

//process merge_bam {
//    tag "$bam.simpleName"
//    publishDir params.output
//    
//    input:
//      path '*.rg.bam' from bams
//
//    output:
//      path "${bam}.merged.bam"
//
//    script:
//    """
//    "echo *.rg.bam"
//   """
//}
