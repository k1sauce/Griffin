include { GRIFFIN_COVERAGE } from '../../modules/local/GRIFFIN_COVERAGE.nf'


workflow GRIFFIN_NUCLEOSOME_PROFILING {
    
    take:
    input_ch         // tuple val(meta), path(bam), path(bai), path(gc_bias_txt)
    reference_genome // ./Ref/hg38.fa
    chrom_sizes_path // ./Ref/hg38.standard.chrom.sizes
    site_file        // /Users/kyle/Projects/Griffin/demo/griffin_nucleosome_profiling_demo_files/sites/CTCF.hg38.1000.txt
    chroms           // chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22
    norm_window      // -5000 5000
    size_range       // 100 200
    map_quality      // 20

    mappability_bw   // ./Ref/k100_minus_exclusion_lists.mappable_regions.hg38.bw
    save_window      // -1000 1000
    center_window    // -30 30
    fft_window       // -960 960
    fft_index       // 10
    smoothing_length // 165
    exclude_path     // ./griffin/excluded_regions.bed
    step            // 15
    cna_normalization_flag // False 
    individual_flag // False
    smoothing_flag // True
    exclude_outliers_flag // True
    exclude_zero_mappability_flag // True
    number_of_sties // none
    site_name       // CTCF_demo



    main:
    GRIFFIN_COVERAGE(
        input_ch,
        reference_genome,
        chrom_sizes_path,
        site_file,
        chroms,
        norm_window,
        size_range,
        map_quality
    )

    GRIFFIN_MERGE_SITES(
        GRIFFIN_COVERAGE.out.gc_bw,
        mappability_bw,
        chrom_sizes_path,
        site_file,
        chroms,
        norm_window,
        save_window,
        center_window,
        fft_window,
        fft_index,
        smoothing_length,
        exclude_path,
        step,
        cna_normalization_flag,
        individual_flag,
        smoothing_flag,
        exclude_outliers_flag,
        exclude_zero_mappability_flag,
        number_of_sties,
        site_name
    )

    emit:
    gc_corrected_bw = GRIFFIN_COVERAGE.out.gc_corrected_bw
    uncorrected_bw = GRIFFIN_COVERAGE.out.uncorrected_bw
    
}