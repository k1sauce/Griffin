process GRIFFIN_COVERAGE {
    tag "$meta.id"
    
    input:
    tuple val(meta), path(bam), path(bai), path(gc_bias_txt)
    path(reference_genome)
    path(chrom_sizes_path)
    path(site_file) // /Users/kyle/Projects/Griffin/demo/griffin_nucleosome_profiling_demo_files/sites/CTCF.hg38.1000.txt
    val(chroms)
    val(norm_window)
    val(size_range)
    val(map_quality)

    output:
    tuple val(meta), path("*.GC_corrected.bw"), path("*.uncorrected.bw"), emit: gc_bw
    tuple val(meta), path("*.GC_corrected.bw"), emit: gc_corrected_bw
    tuple val(meta), path("*.uncorrected.bw"), emit: uncorrected_bw

    script:
    """
    griffin_coverage.py \\
        --sample_name ${meta.id} \\
        --bam ${bam} \\
        --index_file_path ${bai} \\
        --GC_bias ${gc_bias_txt} \\
        --reference_genome ${reference_genome} \\
        --chrom_sizes_path ${chrom_sizes_path} \\
        --site_file ${site_file} \\
        --chroms ${chroms.join(' ')} \\
        --norm_window ${norm_window.join(' ')}\\
        --size_range ${size_range.join(' ')}\\
        --map_quality ${map_quality}
    """

    stub:
    """
    touch ${meta.id}.GC_corrected.bw
    touch ${meta.id}.uncorrected.bw
    """
}