process CLIPPER {

    publishDir "${params.outdir}/clipper", mode: 'copy'
    container  'docker://brianyee/clipper:5d865bb'
    label 'process_high'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.bed"), emit: bed

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    clipper --species hg19 --bam ${bam} --outfile ${prefix}.clip.peakClusters.bed
    """
}

process CREATEREADNUM {

    publishDir "${params.outdir}/readnum", mode: 'copy'
    container  'docker://staphb/samtools:1.20'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.txt"), emit: readnum

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    samtools view -cF 4 ${bam} > ${prefix}.readnum.txt
    """
}

process REMOVE_UNMAPPED_READS {

    publishDir "${params.outdir}/mapped", mode: 'copy'
    container  'docker://staphb/samtools:1.20'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.bam"), emit: bam

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    samtools view -b -F 4 ${bam} > ${prefix}_mapped.bam
    """
}

process OVERLAP_PEAKS {

    publishDir "${params.outdir}/overlapPeaks", mode: 'copy'
    container  'docker://brianyee/eclip:0.7.0_perl'

    input:
    tuple val(meta), path(background), path(signal)

    output:
    tuple val(meta), path("*.bed"), emit: bed

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    overlap_peakfi_with_bam_PE.pl ${signal.bam} ${background.bam} ${signal.bed} ${signal.readnum} ${background.readnum} ${prefix}.normed.bed
    """
}

process MERGE_OVERLAPPING_PEAKS {

    publishDir "${params.outdir}/compressedOverlapedPeaks", mode: 'copy'
    container  'docker://brianyee/merge_peaks:0.1.0'

    input:
    tuple val(meta), path(bed)

    output:
    tuple val(meta), path("*.bed"), emit: bed

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    compress_l2foldenrpeakfi_for_replicate_overlapping_bedformat_outputfull.pl ${bed} ${prefix}.normed.compressed.bed
    """
}

process MAKE_INFORMATION_CONTENT_FROM_PEAKS {

    publishDir "${params.outdir}/entropy", mode: 'copy'
    container  'docker://brianyee/merge_peaks:0.1.0'

    input:
    tuple val(meta), path(results)

    output:
    tuple val(meta), path("*.full"), emit: full
    tuple val(meta), path("*.bed"), emit: bed

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    make_informationcontent_from_peaks.pl ${results[2]} ${results[0]} ${results[1]} ${prefix}.entropy.full ${prefix}.entropy.excessreads
    full_to_bed.py --input ${prefix}.entropy.full --output ${prefix}.entropy.bed 
    """
}

process IDR {

    publishDir "${params.outdir}/ids", mode: 'copy'
    container  'docker://brianyee/idr:2.0.2'

    input:
    tuple val(meta), path(bed)

    output:
    tuple val(meta), path("*.idr"), emit: idr

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    idr --samples ${bed} --input-file-type bed --rank 5 --peak-merge-method max --plot true --output-file ${prefix}.idr
    """
}

process CREATE_BIGWIG {

    publishDir "${params.outdir}/bigwig", mode: 'copy'
    container  'docker://brianyee/makebigwigfiles:0.0.3'

    input:
    tuple val(meta), path(bam)
    file genome

    output:
    tuple val(meta), path("*.bw"), emit: bw

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    makebigwigfiles --bw_pos ${prefix}.pos.bw --bw_neg ${prefix}.neg.bw --bam ${bam} --genome ${genome} --direction r
    """
}
