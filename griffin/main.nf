#!/usr/bin/env nextflow

params.input = "$baseDir/input/*.bam"
params.output = "$baseDir/results"

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

process GRIFFIN_GC_BIAS {

    tag "$meta.id"

    input:
    tuple val(meta), path(gc_counts)
    val mappable_name // 'k100_minus_exclusion_lists.mappable_regions.hg38'
    val genome_GC_frequency // './Ref/genome_GC_frequency'
    val size_range_low
    val size_range_high // [15, 500]

    output:
    path "*.GC_bias.txt"

    script:
    """
    griffin_GC_bias.py \\
        --gc_counts ${gc_counts} \\
        --bam_file_name ${meta.bam_file_name} \\
        --mappable_name ${mappable_name} \\
        --genome_GC_frequency ${genome_GC_frequency} \\
        --size_range ${size_range_low} ${size_range_high}
    """
}

workflow GRIFFIN {

    // input_ch = Channel.from([
    //     [
    //         [
    //             id: 'Healthy_demo',
    //             bam_file_name: 'Healthy_demo',
    //         ], 
    //         file('/Users/kyle/Projects/Griffin/demo/bam/Healthy_GSM1833219_downsampled.sorted.mini.bam'),
    //         file('/Users/kyle/Projects/Griffin/demo/bam/Healthy_GSM1833219_downsampled.sorted.mini.bam.bai')
    //     ]
    // ])

    mappable_regions_path = file('/Users/kyle/Projects/Griffin/Ref/k100_minus_exclusion_lists.mappable_regions.hg38.bed')
    ref_seq = file('/Users/kyle/Projects/Griffin/Ref/hg38.fa')
    chrom_sizes = file('/Users/kyle/Projects/Griffin/Ref/hg38.standard.chrom.sizes')
    map_q = 20
    size_range_low = 15
    size_range_high = 500

    // GC bias parameters
    mappable_name = 'k100_minus_exclusion_lists.mappable_regions.hg38'
    genome_GC_frequency_ch = Channel.fromPath(
        '/Users/kyle/Projects/Griffin/Ref/genome_GC_frequency', 
        type : 'dir', 
        maxDepth : 0
    )

    // GRIFFIN_GC_COUNTS(
    //     input_ch,
    //     mappable_regions_path,
    //     ref_seq,
    //     chrom_sizes,
    //     map_q,
    //     size_range_low,
    //     size_range_high
    // )
    
    // GRIFFIN_GC_COUNTS.out.gc_counts.view()

    input_ch = Channel.from([
        [
            [
                id: 'Healthy_demo',
                bam_file_name: 'Healthy_demo',
            ],
            file('/Users/kyle/Projects/Griffin/results/GC_counts/Healthy_demo.GC_counts.txt'),
        ]
    ])

    GRIFFIN_GC_BIAS(
        input_ch,
        mappable_name,
        genome_GC_frequency_ch,
        size_range_low,
        size_range_high
    )
}

workflow {
    GRIFFIN()
}