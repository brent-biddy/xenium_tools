
process execute_notebook {

    tag "${sample_name}"
    label 'seurat'

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

process bp_cells_clustering {

    tag "${sample_name}"
    label 'seurat'

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore'}

    stageInMode 'copy'

    input:
    path(notebook_path)
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

process execute_squidpy_marker_workflow {

    tag "${sample_name}"
    label 'squidpy'

    time = { 15.m * (1 + task.attempt) }

    input:
    tuple val(sample_name), path(spatialdata_zarr)
    path(notebook)
    path(marker_yaml, stageAs: "sci_adv_markers.yaml")

    output:
    tuple val(sample_name), path("squidpy_marker_report.html"), emit: html

    publishDir "${params.output_path}/results/${sample_name}",
        pattern: "squidpy_marker_report.html",
        saveAs: { "${sample_name}_squidpy_marker_report.html" },
        mode: 'copy'

    script:
    """
    export XDG_CACHE_HOME="./.xdg_cache_home"
    export XDG_DATA_HOME="./.xdg_data_home"
    quarto render ${notebook} -P xenium_zarr:${spatialdata_zarr} --output squidpy_marker_report.html
    """
    stub:
    """
    touch squidpy_marker_report.html
    """
}

process execute_squidpy_notebook {

    tag "${sample_name}"
    label 'squidpy'

    time = { 15.m * (1 + task.attempt) }

    input:
    tuple val(sample_name), path(spatialdata_zarr)
    path(notebook)

    output:
    tuple val(sample_name), path("squidpy_report.html"), emit: html

    publishDir "${params.output_path}/results/${sample_name}",
        pattern: "squidpy_report.html",
        saveAs: { "${sample_name}_squidpy_report.html" },
        mode: 'copy'

    script:
    """
    export XDG_CACHE_HOME="./.xdg_cache_home"
    export XDG_DATA_HOME="./.xdg_data_home"
    quarto render ${notebook} -P xenium_zarr:${spatialdata_zarr} --output squidpy_report.html
    """
    stub:
    """
    touch squidpy_report.html
    """
}

workflow NOTEBOOK_EXECUTION {
    take:
        seurat_rds_ch    // [sample_name, rds] - used for qc_plots and score_markers
        cluster_rds_ch   // [sample_name, rds] - used for cluster_plots
        xenium_dir_ch    // [sample_name, dir] - used for bp_clustering
        spatialdata_ch   // [sample_name, zarr] - used for squidpy_notebook

    main:
        if (params.run_qc_plots) {
            notebook_file = file("${projectDir}/notebooks/xenium_qc_plots.ipynb")
            execute_notebook(seurat_rds_ch, notebook_file, "xenium_qc_report")
        }

        if (params.run_cluster_plots) {
            cluster_notebook = file("${projectDir}/notebooks/seurat_cluster_plots.ipynb")
            execute_notebook(cluster_rds_ch, cluster_notebook, "seurat_cluster_report")
        }

        if (params.score_markers) {
            marker_notebook = file("${projectDir}/notebooks/marker_scores.ipynb")
            marker_yaml = file(params.marker_yaml)
            seurat_rds_with_yaml = seurat_rds_ch.map{ sample_name, rds -> tuple(sample_name, [marker_yaml, rds]) }
            execute_notebook(seurat_rds_with_yaml, marker_notebook, "celltype_marker_gene_report")
        }

        if (params.bp_clustering) {
            bp_clustering_nb = file("${projectDir}/notebooks/bp_cells_clustering.ipynb")
            bp_cells_clustering(bp_clustering_nb, xenium_dir_ch)
        }

        if (params.run_squidpy_notebook) {
            squidpy_nb = file("${projectDir}/notebooks/squidpy_notebook.qmd")
            execute_squidpy_notebook(spatialdata_ch, squidpy_nb)
        }

        if (params.run_squidpy_marker_workflow) {
            squidpy_marker_nb = file("${projectDir}/notebooks/squidpy_marker_workflow.qmd")
            squidpy_marker_yaml = file(params.squidpy_marker_yaml)
            execute_squidpy_marker_workflow(spatialdata_ch, squidpy_marker_nb, squidpy_marker_yaml)
        }
}
