/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FASTQC                 } from '../modules/nf-core/fastqc/main'
include { CUTADAPT               } from '../modules/nf-core/cutadapt/main'
include { CUTADAPT as CUTADAPT2  } from '../modules/nf-core/cutadapt/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { FASTQ_ALIGN_STAR       } from '../subworkflows/nf-core/fastq_align_star/main'
include { BAM_MARKDUPLICATES_PICARD } from '../subworkflows/nf-core/bam_markduplicates_picard/main' 
include { DEEPTOOLS_BAMCOVERAGE } from '../modules/nf-core/deeptools/bamcoverage/main' 
include { SAMTOOLS_MERGE } from '../modules/nf-core/samtools/merge/main'       
include { REMOVE_UNMAPPED_READS  } from '../modules/local/processes.nf'
include { CLIPPER          } from '../modules/local/processes.nf'
include { CREATEREADNUM          } from '../modules/local/processes.nf'
include { OVERLAP_PEAKS          } from '../modules/local/processes.nf'
include { CREATE_BIGWIG           } from '../modules/local/processes.nf'
include { MERGE_OVERLAPPING_PEAKS          } from '../modules/local/processes.nf'
include { MAKE_INFORMATION_CONTENT_FROM_PEAKS   } from '../modules/local/processes.nf'
include { IDR } from '../modules/local/processes.nf'
include { PARSE_IDR_PEAKS } from '../modules/local/processes.nf'
include { OVERLAP_PEAKS_WITH_IDR     } from '../modules/local/processes.nf'
include { GET_REPRODUCING_PEAKS  } from '../modules/local/processes.nf'
include { paramsSummaryMap       } from 'plugin/nf-validation'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_eclipseq_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

overlapPeaksOutputDir = "overlapPeaks"

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

    CUTADAPT2 (
    	CUTADAPT.out.reads
    )


    FASTQ_ALIGN_STAR(
        CUTADAPT2.out.reads,
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

    REMOVE_UNMAPPED_READS(BAM_MARKDUPLICATES_PICARD.out.bam)

    CREATE_BIGWIG(REMOVE_UNMAPPED_READS.out.bam,file(params.chromosomeSize))

    CLIPPER(REMOVE_UNMAPPED_READS.out.bam,params.species)

    CREATEREADNUM(REMOVE_UNMAPPED_READS.out.bam)
   

    REMOVE_UNMAPPED_READS.out.bam.join(CREATEREADNUM.out.readnum).join(CLIPPER.out.bed)
    .map { result  ->
         tuple(result[0].sample,result[0].replicate,[type:result[0].type,bam:result[1],readnum:result[2],bed:result[3]])
    }
    .groupTuple(
         by: [0, 1],
         sort: { e1, e2 -> e1.type <=> e2.type }
    )
    .map {
       result ->
       [[id:result[0],replicate:result[1]],[result[2][0].bam,result[2][0].readnum,result[2][0].bed],[result[2][1].bam,result[2][1].readnum,result[2][1].bed]]
    }
    .set { ch_bamreadbed }

    OVERLAP_PEAKS(ch_bamreadbed)

    MERGE_OVERLAPPING_PEAKS(OVERLAP_PEAKS.out.bedfull)

    ch_bamreadbed
    .map { result -> 
   	 tuple([id:result[0].id,replicate:result[0].replicate],[result[1][1],result[2][1]])
    }
    .set{ ch_readnum }


    MERGE_OVERLAPPING_PEAKS.out.bedfull.join(ch_readnum)
    .set {ch_compressedFullBedAndReadnum}

    MAKE_INFORMATION_CONTENT_FROM_PEAKS(ch_compressedFullBedAndReadnum)

    MAKE_INFORMATION_CONTENT_FROM_PEAKS.out.bed
    .map { result  ->
         tuple(result[0].id,result[0].replicate,result[1])
    }
    .groupTuple(
         by: [0],
         sort: { e1, e2 -> e1[1] <=> e2[1] }
    )
    .map {
       result ->
       [[id:result[0]],result[2]]
    }
    .set { ch_replicatedBeds }

    IDR(ch_replicatedBeds)

    MAKE_INFORMATION_CONTENT_FROM_PEAKS.out.full
    .map { result -> 
   	 tuple([id:result[0].id],result[1])
    }
    .groupTuple(
         by: [0],
         sort: { e1, e2 -> e2 <=> e1 }
    )
    .set{ ch_entropyFiles }

    IDR.out.join(ch_entropyFiles)
    .set{ ch_idrWithEntropy }

    PARSE_IDR_PEAKS(ch_idrWithEntropy)

    REMOVE_UNMAPPED_READS.out.bam.join(CREATEREADNUM.out.readnum)
    .map { result  ->
         tuple(result[0].sample,result[0].replicate,[type:result[0].type,bam:result[1],readnum:result[2]])
    }
    .groupTuple(
         by: [0, 1],
         sort: { e1, e2 -> e1.type <=> e2.type }
    )
    .map {
       result ->
       tuple(result[0],[result[1],result[2]])
     }
     .set { ch_bamreadidr }

    IDR.out
    .map {
       result -> tuple(result[0].id,result[1]) 
    }
    .set {ch_idrWithKey}

    ch_idrWithKey.cross(ch_bamreadidr)
    .map {
       result -> tuple([id:result[0][0],replicate:result[1][1][0]],[result[1][1][1][0].bam,result[1][1][1][0].readnum],[result[1][1][1][1].bam,result[1][1][1][1].readnum,result[0][1]]) 
    } 
    .set {ch_bamreadidr}

    OVERLAP_PEAKS_WITH_IDR(ch_bamreadidr)
    
    OVERLAP_PEAKS_WITH_IDR.out.bedfull.join(MAKE_INFORMATION_CONTENT_FROM_PEAKS.out.full)
    .map {
       result -> tuple(result[0].id,[replicate:result[0].replicate,full:result[1],entropy:result[2]])
    }
    .groupTuple(
         by: [0],
         sort: { e1, e2 -> e1.replicate <=> e2.replicate }   
    )
    .join(ch_idrWithKey)
    .map { result -> tuple(id:result[0],result[1][0],result[1][1],result[2])}
    .set{ch_fullEntropyIDR}

    GET_REPRODUCING_PEAKS(ch_fullEntropyIDR)

    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
    ch_multiqc_files = ch_multiqc_files.mix(CUTADAPT.out.log.collect{it[1]})
    ch_multiqc_files = ch_multiqc_files.mix(FASTQ_ALIGN_STAR.out.log_final.collect{it[1]})
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
