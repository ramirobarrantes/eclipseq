process CLIPPER {

    publishDir "${params.outdir}/clipper", mode: 'copy'
    container  'docker://brianyee/clipper:5d865bb'
    label 'process_high'

    input:
    tuple val(meta), path(bam)
    val species

    output:
    tuple val(meta), path("*.bed"), emit: bed

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    clipper --species ${species} --bam ${bam} --outfile ${prefix}.clip.peakClusters.bed
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
    tuple val(meta), path("*.bed.full"), emit: bedfull

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def replicate = task.ext.prefix ?: "${meta.replicate}"
    """
    overlap_peakfi_with_bam_PE.pl ${signal[0]} ${background[0]} ${signal[2]} ${signal[1]} ${background[1]} ${prefix}_${replicate}.normed.bed
    """
}

process OVERLAP_PEAKS_WITH_IDR {

    publishDir "${params.outdir}/overlapPeaksWithIDR", mode: 'copy'
    container  'docker://brianyee/eclip:0.7.0_perl'

    input:
    tuple val(meta), path(background), path(signal)

    output:
    tuple val(meta), path("*.bed"), emit: bed
    tuple val(meta), path("*.bed.full"), emit: bedfull

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def replicate = task.ext.prefix ?: "${meta.replicate}"
    """
    overlap_peakfi_with_bam_PE.pl ${signal[0]} ${background[0]} ${signal[2]} ${signal[1]} ${background[1]} ${prefix}_${replicate}.idr.normed.bed
    """
}


process MERGE_OVERLAPPING_PEAKS {

    publishDir "${params.outdir}/compressedOverlapedPeaks", mode: 'copy'
    container  'docker://brianyee/merge_peaks:0.1.0'

    input:
    tuple val(meta), path(bed)

    output:
    tuple val(meta), path("*.bed"), emit: bed
    tuple val(meta), path("*.bed.full"), emit: bedfull

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def replicate = task.ext.prefix ?: "${meta.replicate}"
    """
    compress_l2foldenrpeakfi_for_replicate_overlapping_bedformat_outputfull.pl ${bed} ${prefix}_${replicate}.normed.compressed.bed ${prefix}_${replicate}.normed.compressed.bed.full
    """
}

process MAKE_INFORMATION_CONTENT_FROM_PEAKS {

    publishDir "${params.outdir}/entropy", mode: 'copy'
    container  'docker://brianyee/merge_peaks:0.1.0'

    input:
    tuple val(meta), path(compressed_bed), path(readnum)

    output:
    tuple val(meta), path("*.full"), emit: full
    tuple val(meta), path("*.excessreads"), emit: excessreads
    tuple val(meta), path("*.bed"), emit: bed

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def replicate = task.ext.prefix ?: "${meta.replicate}"
    """
    make_informationcontent_from_peaks.pl ${compressed_bed} ${readnum[1]} ${readnum[0]} ${prefix}_${replicate}.entropy.full ${prefix}_${replicate}.entropy.excessreads   
    full_to_bed.py --input ${prefix}_${replicate}.entropy.full --output ${prefix}_${replicate}.entropy.bed
    """
}

process IDR {

    publishDir "${params.outdir}/idr", mode: 'copy'
    container  'docker://brianyee/idr:2.0.2'

    input:
    tuple val(meta), path(bed)

    output:
    tuple val(meta), path("*.idr"), emit: idr

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    idr --samples ${bed[0]} ${bed[1]} --input-file-type bed --rank 5 --peak-merge-method max --plot  --output-file ${prefix}.idr
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


process PARSE_IDR_PEAKS {

    publishDir "${params.outdir}/parseIDRPeaks", mode: 'copy'
    container  'docker://brianyee/merge_peaks:0.1.0'

    input:
    tuple val(meta), path(idr), path(entropy)

    output:
    tuple val(meta), path("*.txt"), emit: summaryStats

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    parse_idr_peaks.pl ${idr} ${entropy[0]} ${entropy[1]} ${prefix}.summaryStats.txt
    """
}


process GET_REPRODUCING_PEAKS {

    publishDir "${params.outdir}/reproducingPeaks", mode: 'copy'
    container  'docker://brianyee/merge_peaks:0.1.0'

    input:
    tuple val(meta), val(replicate1), val(replicate2), path(idr) 

    output:
    tuple val(meta), path("*.full"), emit: full
    tuple val(meta), path("*.bed"), emit: bed

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    get_reproducing_peaks.pl \
       ${replicate1.full} \
       ${replicate2.full} \
       ${prefix}.${replicate1.replicate}.final.full \
       ${prefix}.${replicate2.replicate}.final.full \
       ${prefix}.final.bed \
       ${prefix}.final.custom.bed \
       ${replicate1.entropy} \
       ${replicate2.entropy} ${idr}
    """ 
}

