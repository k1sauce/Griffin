process GRIFFIN_GC_COUNTS {

    tag "$meta.id"

    input:
    tuple val(meta), path(bam), path(bai)
    val mappable_regions_path
    val ref_seq
    val chrom_sizes
    val map_q
    val size_range_low
    val size_range_high

    output:
    tuple val(meta), path("*.GC_counts.txt"), emit: gc_counts

    script:
    """
    griffin_GC_counts.py \\
        --bam_file ${bam} \\
        --bam_file_name ${meta.bam_file_name} \\
        --index_file_path ${bai} \\
        --mappable_regions_path ${mappable_regions_path} \\
        --ref_seq ${ref_seq} \\
        --chrom_sizes ${chrom_sizes} \\
        --map_q ${map_q} \\
        --size_range ${size_range_low} ${size_range_high} \\
        --CPU 1
    """
}