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
  --exclude_path ./griffin/excluded_regions.tsv \
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
    val identifier

    output:
    val identifier

    script:
    """
    
    """

    stub:
    """
    
    """
}