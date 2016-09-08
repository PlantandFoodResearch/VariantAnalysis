#!/usr/bin/env nextflow

log.info ""
log.info '--------------------------------------------------'
log.info 'NEXTFLOW Variant Calling PAIRED END'
log.info '--------------------------------------------------'

if(params.help){
    log.info ''
    log.info 'Usage: '
    log.info 'nextflow run --fastq_folder FASTQ/ --genome genome.fa --species species --genus genus'
    log.info ''
    log.info 'Mandatory arguments:'
    log.info '    --fastq_folder   FOLDER                  Folder containing FASTQ paired files'
    log.info '    --genome FILE                            File containing the reference genome'
    log.info '    --species STRING                         The species you are investigating e.g. chinensis'
    log.info '    --genus String                           The genus (first name) of the species e.g. Actinidia'
    log.info 'Options:'
    log.info ''
    exit 1
}


/*
 * Files
 */

/*
 * The reference genome file
 * Config in nextflow.config:
 */

label_variant_calling_gatk = 'variant_calling_gatk'


genome_file = file(params.genome)


/*
 * Creates the `read_pairs` channel that emits for each read-pair a tuple containing
 * three elements: the pair ID, the first read-pair file and the second read-pair file
 * This is the start of the pipeline and needs to be configured in nextflow.config (params.read).
 *
 * The raw read files
 * Config in nextflow.config:
 */

Channel
    .fromPath("${params.output_dir}/240.add_read_group_id/*.bam")
    .set{ aligned_bam }

Channel
    .fromPath("${params.output_dir}/240.add_read_group_id/*.bai")
    .set{ aligned_bai }


label_variant_calling_gatk = 'variant_calling_gatk'


process index_faidx_unzip{
    input:
    file genome from genome_file

    output:
    file 'kiwitest.fasta' into genome_raw

    """
    gunzip -c ${genome} > kiwitest.fasta
    """

}


process gatk_haplotype_calling{

    module = params.variant_calling_gatk_module

    publishDir "${params.publish_dir}/260.${tag}"

    input:
    file genome from genome_raw
    file bwa_mdup_rg_bam  from aligned_bam
    file bwa_mdup_rg_bai  from aligned_bai

    output:
    file "${bwa_mdup_rg_bam}.vcf" into gatk_gvcf

    script:
    """
    samtools faidx ${genome}

    tt=${genome}

    java -jar \$PICARD CreateSequenceDictionary \
            R=${genome} O=\${tt%.fasta}.dict

    java -Djava.io.tmpdir=/workspace/cfphxd/tmp -jar /software/bioinformatics/gatk-1.0/GenomeAnalysisTK.jar  \
      -T HaplotypeCaller \
      -R ${genome} \
      -I ${bwa_mdup_rg_bam} \
      --genotyping_mode DISCOVERY \
      -ploidy ${params.ploidy} \
      -o ${bwa_mdup_rg_bam}.vcf

      ln -s ${params.publish_dir}/260.${label_variant_calling_gatk} ${params.ouput_dir}
    """
}


gvcf_list = file("VCFlist.list")

process gatk_joint_calling{

    module = params.joint_calling_module

    input:
    file genome from genome_file
    file 'vcf_*' from gatk_gvcf.toList()
    file gvcf_list from gvcf_list

    output:
    file 'gvcf_joint.g.vcf' into gvcf_joint

    publishDir "${params.publish_dir}/270.joint_calling"

    script:
    """
    samtools faidx ${genome}

    tt=${genome}

    java -jar \$PICARD CreateSequenceDictionary \
        R=${genome} O=\${tt%.fasta}.dict

    cat ${gvcf_list} > list.list
    ls vcf_* >> list.list

    java -Djava.io.tmpdir=/workspace/cfphxd/tmp -jar /software/bioinformatics/gatk-1.0/GenomeAnalysisTK.jar \
        -T GenotypeGVCFs \
        -nt ${params.joint_calling_nt} \
        -R ${genome} \
        -V list.list \
        -o gvcf_joint.g.vcf

        ln -s ${params.publish_dir}/270.joint_calling ${params.ouput_dir}
    """
}


