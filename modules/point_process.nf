process bp_cells_clustering {

    tag "${sample_name}"

    // time { 2.h * (1 + task.attempt)}
    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore'}

    stageInMode 'copy'

    input:
    path (notebook_path)
    tuple val(sample_name), path(xenium_output)

    output:
    tuple val(sample_name), path("jupyter_notebook.html"), emit: html
    tuple val(sample_name), path("clustering_results.csv"), emit: cluster_csv, optional: true

    publishDir "${params.output_path}/results/${sample_name}/bp_cells/", pattern: "jupyter_notebook.html", saveAs: { "${sample_name}_bp_cells_clustering_report.html" }, mode: 'copy'
    publishDir "${params.output_path}/results/${sample_name}/bp_cells/", pattern: "clustering_results.csv", saveAs: { "${sample_name}_bp_cells_clustering_results.csv" }, mode: 'copy'

    script:
    """
    jupyter nbconvert --execute --allow-errors --output jupyter_notebook --to html ${notebook_path}
    """
    stub:
    """
    touch jupyter_notebook.html
    """
}

workflow TRANSCRIPTS {
    take:
        sample_info

    main:

    sample_info.view()
}