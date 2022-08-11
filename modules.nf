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
    publishDir params.results

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
    publishDir params.results 

    input:
      path sam

    output:
      path "${sam.baseName}.sorted.bam"

    script:
    """
    picard SortSam -I ${sam} -O ${sam.baseName}.sorted.bam -SO coordinate --VALIDATION_STRINGENCY LENIENT --CREATE_INDEX TRUE
    """
}

process fix_mate {
    tag "$bam.simpleName"
    
    input:
      path bam

    output:
      path "${bam.simpleName}.fixmate.bam"

    script:
    """
    picard FixMateInformation -I ${bam} -O ${bam.simpleName}.fixmate.bam --VALIDATION_STRINGENCY LENIENT --CREATE_INDEX TRUE
    """
}
