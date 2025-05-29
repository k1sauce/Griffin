process GRIFFIN_PLOT {
    tag "$meta.id"

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