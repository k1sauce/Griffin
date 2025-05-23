// ./griffin/bin/griffin_plot.py \
//     --sample Healthy_demo \
//     --uncorrected_coverage_tsv Healthy_demo.uncorrected.coverage.tsv \
//     --gc_corrected_coverage_tsv Healthy_demo.GC_corrected.coverage.tsv \
//     --save_window -1000 1000 \
//     --step 15 

process GRIFFIN_PLOT {
    input:
    tuple val(meta), path(uncorrected_coverage_tsv), path(gc_corrected_coverage_tsv)
    val(save_window)
    val(step)

    output:
    path("*.pdf"), emit: pdf

    script:
    """
    griffin_plot.py \\
        --sample ${meta.id} \\
        --uncorrected_coverage_tsv ${uncorrected_coverage_tsv} \\
        --gc_corrected_coverage_tsv ${gc_corrected_coverage_tsv} \\
        --save_window ${save_window.join(' ')} \\
        --step ${step} 
    """

    stub:
    """
    touch CTCF_demo.pdf
    """
}