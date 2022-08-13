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
    bwa mem \
        -M \
        ${genome} \
        ${reads[0]} \
        ${reads[1]} \
        > ${sampleName}.sam
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
    gatk SortSam \
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
    gatk FixMateInformation \
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
    gatk MarkDuplicates \
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
    gatk AddOrReplaceReadGroups \
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

process merge_bam {
    publishDir params.output, mode: 'copy'
    
    input:
        path params.output
        path add_read_groups

    output:
        path "merged.bam"
    
    script:
    """
    samtools merge \
        merged.bam \
        $params.output/*.rg.bam
    """
}

process haplotype_caller {
    tag "$merged_bam.baseName"
    publishDir params.output, mode: 'copy'
    
    input:
        path genome
        path genome_fai
        path genome_dic
        path merged_bam
      
    output:
        path "${merged_bam.baseName}.vcf.gz"

    script:
    """
    samtools index ${merged_bam}

    gatk HaplotypeCaller \
      -R ${genome} \
      -I ${merged_bam} \
      -O ${merged_bam.baseName}.vcf.gz \
      -ERC GVCF
    """
}
