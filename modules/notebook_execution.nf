
process execute_notebook {

    tag "${sample_name}"

    time = { 15.m * (1 + task.attempt) }

    input:
    tuple val(sample_name), path(input_files)
    path(notebook)
    val(output_name)

    output:
    tuple val(sample_name), path("jupyter_notebook.html"), emit: html

    publishDir "${params.output_path}/results/${sample_name}",
        pattern: "jupyter_notebook.html",
        saveAs: { "${sample_name}_${output_name}.html" },
        mode: 'copy'

    script:
    """
    jupyter nbconvert --execute --allow-errors --output jupyter_notebook --to html ${notebook}
    """
    stub:
    """
    touch jupyter_notebook.html
    """
}
