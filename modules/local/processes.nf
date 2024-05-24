process CLIPPER {

    publishDir "${params.outdir}/clipper", mode: 'copy'
    container  'docker://brianyee/clipper:5d865bb'

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
