/*
./griffin/bin/griffin_merge_sites.py \
  --sample_name Healthy_demo \
  --uncorrected_bw_path ./tmp/Healthy_demo/tmp_bigWig/Healthy_demo.uncorrected.bw \
  --GC_corrected_bw_path ./tmp/Healthy_demo/tmp_bigWig/Healthy_demo.GC_corrected.bw \
  --mappability_bw ./Ref/k100.Umap.MultiTrackMappability.bw \
  --chrom_sizes_path ./Ref/hg38.standard.chrom.sizes \
  --site_file ./demo/griffin_nucleosome_profiling_demo_files/sites/CTCF.hg38.1000.txt \
  --chroms chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22\
  --norm_window -5000 5000 \
  --save_window -1000 1000 \
  --center_window -30 30 \
  --fft_window -960 960 \
  --fft_index 10 \
  --smoothing_length 165 \
  --exclude_path ./griffin/excluded_regions.bed \
  --step 15 \
  --CNA_normalization False \
  --individual False \
  --smoothing True \
  --exclude_outliers True \
  --exclude_zero_mappability True \
  --number_of_sites none \
  --site_name CTCF_demo
*/

process GRIFFIN_MERGE_SITES {
    tag "$meta.id"

    input:
    tuple val(meta), path(corrected_bw), path(uncorrected_bw)
    path(mappability_bw)
    path(chrom_sizes_path)
    path(site_file) 
    val(chroms)
    val(norm_window)
    val(save_window)
    val(center_window)
    val(fft_window)
    val(fft_index)
    val(smoothing_length)
    path(exclude_path)
    val(step)
    val(cna_normalization_flag)
    val(individual_flag)
    val(smoothing_flag)
    val(exclude_outliers_flag)
    val(exclude_zero_mappability_flag)
    val(number_of_sties)
    val(site_name)

    output:
    tuple val(meta), path("${meta.id}.GC_corrected.coverage.tsv"), path("${meta.id}.uncorrected.coverage.tsv"), emit: merge_sites

    script:
    """
    griffin_merge_sites.py \\
        --sample_name ${meta.id} \\
        --uncorrected_bw_path ${uncorrected_bw} \\
        --GC_corrected_bw_path ${corrected_bw} \\
        --mappability_bw ${mappability_bw} \\
        --chrom_sizes_path ${chrom_sizes_path} \\
        --site_file ${site_file} \\
        --chroms ${chroms.join(' ')} \\
        --norm_window ${norm_window.join(' ')} \\
        --save_window ${save_window.join(' ')} \\
        --center_window ${center_window.join(' ')} \\
        --fft_window ${fft_window.join(' ')} \\
        --fft_index ${fft_index} \\
        --smoothing_length ${smoothing_length} \\
        --exclude_path ${exclude_path} \\
        --step ${step} \\
        --CNA_normalization ${cna_normalization_flag ? 'True' : 'False'} \\
        --individual ${individual_flag ? 'True' : 'False'} \\
        --smoothing ${smoothing_flag ? 'True' : 'False'} \\
        --exclude_outliers ${exclude_outliers_flag ? 'True' : 'False'} \\
        --exclude_zero_mappability ${exclude_zero_mappability_flag ? 'True' : 'False'} \\
        --number_of_sites ${number_of_sties} \\
        --site_name ${site_name}

    """

    stub:
    """
    touch ${meta.id}.Healthy_demo.GC_corrected.coverage.tsv
    touch ${meta.id}.Healthy_demo.uncorrected.coverage.tsv
    """
}