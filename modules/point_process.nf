process point_process {

    tag "${sample_name}"

    // time { 2.h * (1 + task.attempt)}
    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore'}

    container "babiddy755/xenium_tools_pasta"

    stageInMode 'copy'

    input:
    path (notebook_path)
    path (notebook_path_2)
    path (script_path)
    tuple val(sample_name), path(xenium_output)

    output:
    tuple val(sample_name), path("jupyter_notebook.html"), emit: html
    tuple val(sample_name), path("clustering_results.csv"), emit: cluster_csv, optional: true

    publishDir "${params.output_path}/results/${sample_name}/bp_cells/", pattern: "jupyter_notebook.html", saveAs: { "${sample_name}_bp_cells_clustering_report.html" }, mode: 'copy'
    publishDir "${params.output_path}/results/${sample_name}/bp_cells/", pattern: "clustering_results.csv", saveAs: { "${sample_name}_bp_cells_clustering_results.csv" }, mode: 'copy'

    script:
    """
    jupyter nbconvert --execute --allow-errors --output jupyter_notebook --to html ${notebook_path}
    R -e "rmarkdown::render('${notebook_path_2}')"
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

    point_process_notebook = file("${projectDir}/notebooks/setup.ipynb")
    point_process_notebook_2 = file("${projectDir}/notebooks/theory_point.Rmd")
    script_file = file("${projectDir}/notebooks/utils.R")
    point_process(point_process_notebook, point_process_notebook_2, script_file, sample_info)
}