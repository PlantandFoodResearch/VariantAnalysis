#!/usr/bin/env nextflow

log.info '--------------------------------------------------'
log.info 'NEXTFLOW Variant Calling PAIRED END'
log.info '--------------------------------------------------'

if(params.help){
    log.info ''
    log.info 'Usage: '
    log.info 'TEST: nextflow run hdzierz/VariantAnalysis/align.nf'
    log.info 'PROD: nextflow run hdzierz/VariantAnalysis/align.nf --input_dir FASTQ/ --genome genome.fa --species species --genus genus'
    log.info ''
    log.info 'Mandatory arguments:'
    log.info '    --fastq_folder   FOLDER                  Folder containing FASTQ paired files'
    log.info '    --genome FILE                            File containing the reference genome'
    log.info '    --species STRING                         The species you are investigating e.g. chinensis'
    log.info '    --genus String                           The genus (first name) of the species e.g. Actinidia'
    log.info 'Options:
    log.info '    --design FILE                            experimental design description file'
    log.info ''
    exit 1
}

