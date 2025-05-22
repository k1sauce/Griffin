#!/usr/bin/env nextflow
include { GRIFFIN_GC_COUNTS } from '../../modules/local/GRIFFIN_GC_COUNTS.nf'
include { GRIFFIN_GC_BIAS } from '../../modules/local/GRIFFIN_GC_BIAS.nf'

workflow GRIFFIN_GC_AND_MAPPABILITY_CORRECTION {
    take:
    input_ch // tuple val(meta), path(bam), path(bai)
    mappable_regions_path // '/Users/kyle/Projects/Griffin/Ref/k100_minus_exclusion_lists.mappable_regions.hg38.bed'
    ref_seq // '/Users/kyle/Projects/Griffin/Ref/hg38.fa'
    chrom_sizes // '/Users/kyle/Projects/Griffin/Ref/hg38.standard.chrom.sizes'
    map_q // 20 
    size_range_low // 15
    size_range_high // 500
    mappable_name // 'k100_minus_exclusion_lists.mappable_regions.hg38'
    genome_GC_frequency_ch // Channel.fromPath('/Users/kyle/Projects/Griffin/Ref/genome_GC_frequency', type : 'dir', maxDepth : 0)

    main:

    GRIFFIN_GC_COUNTS(
        input_ch,
        mappable_regions_path,
        ref_seq,
        chrom_sizes,
        map_q,
        size_range_low,
        size_range_high
    )

    GRIFFIN_GC_BIAS(
        GRIFFIN_GC_COUNTS.out.gc_counts,
        mappable_name,
        genome_GC_frequency_ch,
        size_range_low,
        size_range_high
    )

    emit:
    gc_bias_txt = GRIFFIN_GC_BIAS.out.gc_bias_txt
    gc_bias_pdf = GRIFFIN_GC_BIAS.out.gc_bias_pdf
}