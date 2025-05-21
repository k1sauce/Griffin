#!/usr/bin/env nextflow
include { GRIFFIN_GC_COUNTS } from './modules/local/GRIFFIN_GC_COUNTS.nf'
include { GRIFFIN_GC_BIAS } from './modules/local/GRIFFIN_GC_BIAS.nf'

workflow GRIFFIN_GC_AND_MAPPABILITY_CORRECTION {

    main:
    // This is for running starting from GRIFFIN_GC_COUNTS but 
    // its a bit slow
    //
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

    // Additional GC bias parameters
    mappable_name = 'k100_minus_exclusion_lists.mappable_regions.hg38'
    genome_GC_frequency_ch = Channel.fromPath(
        '/Users/kyle/Projects/Griffin/Ref/genome_GC_frequency', 
        type : 'dir', 
        maxDepth : 0
    )

    GRIFFIN_GC_BIAS(
        input_ch,
        mappable_name,
        genome_GC_frequency_ch,
        size_range_low,
        size_range_high
    )

    emit:
    gc_bias_txt = GRIFFIN_GC_BIAS.out.gc_bias_txt
    gc_bias_pdf = GRIFFIN_GC_BIAS.out.gc_bias_pdf
}


workflow GRIFFIN {
    GRIFFIN_GC_AND_MAPPABILITY_CORRECTION()

}

workflow {
    GRIFFIN()
}