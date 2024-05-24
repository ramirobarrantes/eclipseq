/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FASTQC                 } from '../modules/nf-core/fastqc/main'
include { CUTADAPT               } from '../modules/nf-core/cutadapt/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { FASTQ_ALIGN_STAR       } from '../subworkflows/nf-core/fastq_align_star/main'
include { BAM_MARKDUPLICATES_PICARD } from '../subworkflows/nf-core/bam_markduplicates_picard/main' 
include { DEEPTOOLS_BAMCOVERAGE } from '../modules/nf-core/deeptools/bamcoverage/main' 
include { SAMTOOLS_MERGE } from '../modules/nf-core/samtools/merge/main'       
include { CLIPPER          } from '../modules/local/processes.nf'
include { CREATEREADNUM          } from '../modules/local/processes.nf'
include { paramsSummaryMap       } from 'plugin/nf-validation'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_eclipseq_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow ECLIPSEQ {

    take:
    ch_samplesheet // channel: samplesheet read in from --input

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    //
    // MODULE: Run FastQC
    //
    FASTQC (
        ch_samplesheet
    )

    CUTADAPT (
    	ch_samplesheet
    )

    FASTQ_ALIGN_STAR(
        CUTADAPT.out.reads,
        [[:],file(params.star_index)],
        [[:],file(params.star_gtf)],
        'TRUE',
        'Illumina',
        '',
        [[:],file(params.star_fasta)],
        ''
    )

    BAM_MARKDUPLICATES_PICARD(
	FASTQ_ALIGN_STAR.out.bam,
	[[:],file(params.star_fasta)],
        [[:],file(params.star_fasta_fai)])

    //Merge technical replicates of the same type
    BAM_MARKDUPLICATES_PICARD.out.bam.map {
        meta, bam ->
        fmeta = meta
 	fmeta.newid = fmeta.sample
        [fmeta, bam] }.view()
//        [ fmeta, bam ] }.view()
/*        .groupTuple(by: [0])
        .map { it ->  [ it[0], it[1].flatten() ] }
        .set { ch_merged_bam }
*/   
//     SAMTOOLS_MERGE(ch_merged_bam)

    //CLIPPER(BAM_MARKDUPLICATES_PICARD.out.bam)

//    DEEPTOOLS_BAMCOVERAGE(
//	BAM_MARKDUPLICATES_PICARD.out.bam.join(BAM_MARKDUPLICATES_PICARD.out.bai,by: [0]),
//	file(params.star_fasta),
//        file(params.star_fasta_fai)
//	)

    //CREATEREADNUM(BAM_MARKDUPLICATES_PICARD.out.bam)
   
    //BAM_MARKDUPLICATES_PICARD.out.bam.mix(CREATEREADNUM.out.readnum).map {
    //    meta, result ->
    //    fmeta = meta
    //    fmeta.id = fmeta.id
   //     [ fmeta, result ] }
   //     .groupTuple(by: [0])
    //    .map { it ->  [ it[0], it[1].flatten() ] }
    //    .set { ch_ }

    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
    ch_multiqc_files = ch_multiqc_files.mix(CUTADAPT.out.log.collect{it[1]})
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_pipeline_software_mqc_versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))

    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
