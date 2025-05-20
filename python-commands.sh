## Griffing GC and mappability correction scripts
./scripts/griffin_GC_counts.py \
  --bam_file ./demo/bam/Healthy_GSM1833219_downsampled.sorted.mini.bam \
  --bam_file_name Healthy_demo \
  --mappable_regions_path ./Ref/k100_minus_exclusion_lists.mappable_regions.hg38.bed \
  --ref_seq ./Ref/hg38.fa \
  --chrom_sizes ./Ref/hg38.standard.chrom.sizes \
  --out_dir results \
  --map_q 20 \
  --size_range 15 500 \
  --CPU 1

./scripts/griffin_GC_bias.py \
  --bam_file_name Healthy_demo \
  --mappable_name k100_minus_exclusion_lists.mappable_regions.hg38 \
  --genome_GC_frequency ./Ref/genome_GC_frequency \
  --out_dir results \
  --size_range 15 500

## Griffin nucleosome profiling scripts
./scripts/griffin_coverage.py \
  --sample_name Healthy_demo \
  --bam demo/bam/Healthy_GSM1833219_downsampled.sorted.mini.bam \
  --GC_bias demo/griffin_GC_correction_demo_files/expected_results/Healthy_demo.GC_bias.txt \
  --mappability_bias none \
  --mappability_correction False \
  --tmp_dir tmp \
  --reference_genome ./Ref/hg38.fa \
  --mappability_bw ./Ref/k100.Umap.MultiTrackMappability.bw \
  --chrom_sizes_path ./Ref/hg38.standard.chrom.sizes \
  --sites_yaml ./sites.yaml \
  --griffin_scripts ./scripts \
  --chrom_column Chrom \
  --position_column position \
  --strand_column Strand \
  --chroms chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22\
  --norm_window -5000 5000 \
  --size_range 100 200 \
  --map_quality 20 \
  --number_of_sites none \
  --sort_by none \
  --ascending none \
  --CPU 1 

./scripts/griffin_merge_sites.py \
  --sample_name Healthy_demo \
  --uncorrected_bw_path ./tmp/Healthy_demo/tmp_bigWig/Healthy_demo.uncorrected.bw \
  --GC_corrected_bw_path ./tmp/Healthy_demo/tmp_bigWig/Healthy_demo.GC_corrected.bw \
  --GC_map_corrected_bw_path none \
  --mappability_correction False \
  --tmp_dir ./tmp \
  --results_dir ./results \
  --mappability_bw ./Ref/k100.Umap.MultiTrackMappability.bw \
  --chrom_sizes_path ./Ref/hg38.standard.chrom.sizes \
  --sites_yaml ./sites.yaml \
  --griffin_scripts ./scripts \
  --chrom_column Chrom \
  --position_column position \
  --strand_column Strand \
  --chroms chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22\
  --norm_window -5000 5000 \
  --save_window -1000 1000 \
  --center_window -30 30 \
  --fft_window -960 960 \
  --fft_index 10 \
  --smoothing_length 165 \
  --exclude_paths ./Ref/encode_unified_GRCh38_exclusion_list.bed ./Ref/hg38_centromeres.bed ./Ref/hg38_gaps.bed ./Ref/hg38_fix_patches.bed ./Ref/hg38_alternative_haplotypes.bed\
  --step 15 \
  --CNA_normalization False \
  --individual False \
  --smoothing True \
  --exclude_outliers True \
  --exclude_zero_mappability True \
  --number_of_sites none \
  --sort_by none \
  --ascending none \
  --CPU 1 

./scripts/griffin_plot.py \
    --in_dir results \
    --samples_yaml ./samples.GC.yaml \
    --mappability_correction False \
    --save_window -1000 1000 \
    --step 15 \
    --individual False \
    --out_dir results
    
# find tmp -type d -empty -delete