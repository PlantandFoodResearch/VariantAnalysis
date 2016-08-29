#!/usr/bin/env nextflow

if(!file(params.design).exists()){
    println("COPYING design file to local directory")
    println("cp $baseDir/design.config ${params.design}")
    ['cp', "$baseDir/design.config", params.design].execute().waitFor()
}
else{
    println("Design file exists")
}


f = file("${params.design}")

reader = f.newReader()
f.withReader{
    String line
    while( line = reader.readLine()){
        def row = line.tokenize(',')
        println("Parsing row: $row")
        if(row[0] != 'label'){
            println("Linking files")
            fn = "${params.input_dir}/${row[1]}"
            sl = "${params.data_dir}/${row[0]}_${row[2]}_${row[3]}.fq.gz"
            println("ln -s $fn $sl")
            ['ln', '-s', fn, sl].execute().waitFor()
        }
    }
}


/*
 * Files
 */

/*
 * The reference genome file
 * Config in nextflow.config:
 */

label_align = 'align'
label_mark_dup = 'mark_dup'
label_add_read_group_id = 'add_read_group_id'
label_variant_calling = 'variant_calling'

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
    .fromFilePairs( params.reads )
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
    .set { read_pairs }


/*
 * Indexing target genome. 
 * Config in nextflow.config:
 *  tool: bwa
 *  input: File genome_file
 *  output: Channel genome_index
 *  params.build_index: true
 */

process build_index {
    tag 'build_index'
    module = params.build_index_module

    //beforeScript "source init.sh"

    input:
    file genome_file from genome_file
      
    output:
    file 'genome.index*' into genome_index
       
    when:
    params.build_index == true

    script:
    """
        bwa index -a bwtsw ${genome_file} -p genome.index
    """
}


/*
 * Aligning raw reads against reference genome
 * Config in nextflow.config:
 *  tool: bwa
 *  input: File genoem_file, Channel read_pairs
 *  output: file bam
 */

process align{
    tag label_align
    module = params.align_module

    publishDir "${params.output_dir}/220.${tag}", mode: 'copy', overwrite: true

    input:
    file genome_file from genome_file
    file genome_index from genome_index.first()
    set pair_id, file(reads) from read_pairs

    output:
    set pair_id, file("${tag}_${pair_id}.bam") into bam_align

    when:
    params.align == true

    script:
    """
        bwa mem -t 8 -M genome.index ${reads} | samtools view -u - | samtools sort -O bam -T genome_file -o ${label_align}_${pair_id}.bam
    """
}

/*
 
  Aligning raw reads against reference genome
  Config in nextflow.config
   tool: samtools
   input: Channel bam
   output: File mbam 
 */

/*
process align_merge{
    tag 'align_merge'
    module = params.align_merge_module 

    input:
    file bam_file from bam.toList()

    output:
    file 'merged.bam' into mbam

    when:
    params.align_merge == true

    script:
    """
        samtools merge merged.bam $bam_file
    """
}
*/


/*
 * Mark duplicates
 * Tool: PICARD 
 * Input: Channel mbam
 * Output 
 */

process mark_dup{
    tag label_mark_dup
    module = params.mark_dup_module

    input:
    set pair_id, file(bwa_mem) from bam_align

    output:
    set pair_id, file("${tag}_${pair_id}.bam") into bwa_mdup_bam

    publishDir "${params.output_dir}/230.${tag}", mode: 'copy', overwrite: true

    when:
    params.mark_dup == true

    script:
    """
       java -jar -Xmx32G \$PICARD MarkDuplicates \
       INPUT=$bwa_mem \
       OUTPUT=${label_mark_dup}_${pair_id}.bam \
       AS=true \
       MAX_RECORDS_IN_RAM=50000000 \
       MAX_FILE_HANDLES=1000 \
       M=${label_mark_dup}_${pair_id}.txt \
       OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
       TMP_DIR=tmp
    """
        
}


/*
 * Add Read Group
 * Tool: picard
 * Input: File bwa_mdup_bam
 * Ouput: File bwa_mdup_rg_bam
 */

process add_read_group_id{
    tag label_add_read_group_id

    module = params.add_read_group_id_module

    input:
    set pair_id, file(bwa_mdup_bam) from bwa_mdup_bam

    output:
    set pair_id, file("${tag}_${pair_id}.bam") into bwa_mdup_rg_bam

    publishDir "${params.output_dir}/240.${tag}", mode: 'copy', overwrite: true

    when:
    params.add_read_group_id == true

    script:
    """
        java -jar \$PICARD AddOrReplaceReadGroups \
        I=${bwa_mdup_bam} \
        O=${label_add_read_group_id}_${pair_id}.bam \
        MAX_RECORDS_IN_RAM=${params.add_read_group_id_MAX_RECORDS_IN_RAM} \
        RGID=${pair_id} \
        RGLB=${params.add_read_group_id_RGLB} \
        RGPL=${params.add_read_group_id_RGPL} \
        RGPU=${params.add_read_group_id_RGPU} \
        RGSM=${pair_id}
        ${params.add_read_group_id_free}

        samtools index ${label_add_read_group_id}_${pair_id}.bam
        samtools flagstat ${label_add_read_group_id}_${pair_id}.bam > ${label_add_read_group_id}_${pair_id}.stats
    """
}


/*
 * Variant calling with freebayse
 * Input: Channel bwa_mdup_rg
 * Output: 
 *  File genome
 *  File vcf
 * Tool: freebayse
 */
process variant_calling_freebayse{
    tag label_variant_calling

    module = params.variant_calling_freebayse_module 

    input:
    file genome from genome_file
    set pair_id, file(bwa_mdup_rg_bam) from bwa_mdup_rg_bam

    output:
    set pair_id, file("${tag}_${pair_id}.vcf") into vcf

    publishDir "${params.output_dir}/250.${tag}", mode: 'copy', overwrite: true

    when:
    params.variant_calling_freebayse == true

    script:
    """
       freebayes -f ${genome_file} \
       -b ${bwa_mdup_rg_bam} \
       -v ${label_variant_calling}_${pair_id}.vcf \
       -C 5 \
       --report-genotype-likelihood-max \
       --use-mapping-quality \
       --genotype-qualities -d
    """

}
