process GRIFFIN_GC_BIAS {

    tag "$meta.id"

    input:
    tuple val(meta), path(gc_counts)
    val mappable_name // 'k100_minus_exclusion_lists.mappable_regions.hg38'
    val genome_GC_frequency // './Ref/genome_GC_frequency'
    val size_range_low
    val size_range_high // [15, 500]

    output:
    tuple val(meta), path("*.txt"), emit: gc_bias_txt
    tuple val(meta), path("*.pdf"), emit: gc_bias_pdf

    script:
    """
    griffin_GC_bias.py \\
        --gc_counts ${gc_counts} \\
        --bam_file_name ${meta.bam_file_name} \\
        --mappable_name ${mappable_name} \\
        --genome_GC_frequency ${genome_GC_frequency} \\
        --size_range ${size_range_low} ${size_range_high}
    """

    stub:
    """
    touch ${meta.bam_file_name}.GC_bias.txt
    touch ${meta.bam_file_name}.GC_bias.pdf
    touch ${meta.bam_file_name}.GC_bias.summary.pdf
    touch ${meta.bam_file_name}.GC_bias.key_lengths.pdf
    """
}