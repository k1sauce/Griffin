include { GRIFFIN_GC_AND_MAPPABILITY_CORRECTION } from './subworkflows/local/GRIFFIN_GC_AND_MAPPABILITY_CORRECTION.nf'
include { GRIFFIN_NUCLEOSOME_PROFILING } from './subworkflows/local/GRIFFIN_NUCLEOSOME_PROFILING.nf'

workflow GRIFFIN {

    // Input parameters
    mappable_regions_path = file('/Users/kyle/Projects/Griffin/Ref/k100_minus_exclusion_lists.mappable_regions.hg38.bed')
    ref_seq = file('/Users/kyle/Projects/Griffin/Ref/hg38.fa')
    chrom_sizes = file('/Users/kyle/Projects/Griffin/Ref/hg38.standard.chrom.sizes')
    map_q = 20
    size_range_low = 15
    size_range_high = 500
    mappable_name = 'k100_minus_exclusion_lists.mappable_regions.hg38'
    genome_GC_frequency_dir = '/Users/kyle/Projects/Griffin/Ref/genome_GC_frequency'
    site_file = file('/Users/kyle/Projects/Griffin/demo/griffin_nucleosome_profiling_demo_files/sites/CTCF.hg38.1000.txt')
    chroms = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22']
    normalization_window = [-5000, 5000]
    size_range = [100, 200]
    mappability_bw = file('/Users/kyle/Projects/Griffin/Ref/k100_minus_exclusion_lists.mappable_regions.hg38.bw')
    save_window = [-100, 1000]
    center_window = [-30, 30]
    fft_window = [-960, 960]
    fft_index = 10
    smoothing_length = 165
    exclude_path = '/Users/kyle/Projects/Griffin/griffin/excluded_regions.bed'
    step = 15
    cna_normalization_flag = False
    individual_flag = False
    smoothing_flag = True
    exclude_outliers_flag = True
    exclude_zero_mappability_flag = True
    number_of_sties = 'none' 
    site_name = 'CTCF_demo'

    input_ch = Channel.from([
        [
            [
                id: 'Healthy_demo',
                bam_file_name: 'Healthy_demo',
            ], 
            file('/Users/kyle/Projects/Griffin/demo/bam/Healthy_GSM1833219_downsampled.sorted.mini.bam'),
            file('/Users/kyle/Projects/Griffin/demo/bam/Healthy_GSM1833219_downsampled.sorted.mini.bam.bai')
        ]
    ])

    genome_GC_frequency_ch = Channel.fromPath(
        genome_GC_frequency_dir, 
        type : 'dir', 
        maxDepth : 0
    )

    GRIFFIN_GC_AND_MAPPABILITY_CORRECTION(
        input_ch,
        mappable_regions_path,
        ref_seq,
        chrom_sizes,
        map_q,
        size_range_low,
        size_range_high,
        mappable_name,
        genome_GC_frequency_ch
    )


    // combine input_ch and gc_counts on the ide
    next_ch = input_ch
        .combine(GRIFFIN_GC_AND_MAPPABILITY_CORRECTION.out.gc_bias_txt, by: 0)
        .dump(tag: 'next_ch')

    GRIFFIN_NUCLEOSOME_PROFILING(
        next_ch,
        ref_seq,
        chrom_sizes,
        site_file, 
        chroms,
        normalization_window,
        size_range,
        map_q,
        mappability_bw,   
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
}

workflow {
    GRIFFIN()
}
