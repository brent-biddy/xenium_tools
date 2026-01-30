process cluster_seurat {
    
    tag "${sample_name}"

    input:
    tuple val(sample_name), path(seurat_obj)
    
    output:
    tuple val(sample_name), path("seurat_clusters.RDS"), emit: rds
    tuple val(sample_name), path("seurat_clusters.csv"), emit: csv
    
    publishDir "${params.output_path}/results/${sample_name}", pattern: "seurat_clusters.RDS", saveAs: { "${sample_name}_seurat_clustered.RDS" }, mode: 'copy'
    publishDir "${params.output_path}/results/${sample_name}", pattern: "seurat_clusters.csv", saveAs: { "${sample_name}_seurat_clustered.csv" }, mode: 'copy'


    script:
    """
    cluster_seurat_xenium.R --seurat_object ${seurat_obj} --executor ${task.executor}
    """
    stub:
    """
    touch seurat_clusters.RDS
    touch seurat_clusters.csv
    """
}

process run_notebook {

    tag "${sample_name}"

    input:
    path (notebook_path)
    tuple val(sample_name), path(seurat_rds)

    output:
    tuple val(sample_name), path("jupyter_notebook.html"), emit: html

    publishDir "${params.output_path}/results/${sample_name}", pattern: "jupyter_notebook.html", saveAs: { "${sample_name}_plots.html" }, mode: 'copy'

    script:
    """
    jupyter nbconvert --execute --allow-errors --output jupyter_notebook --to html ${notebook_path}
    """
}

workflow SEURAT {
    take:
        seurat_rds

    main:
        if(params.cluster){
            cluster_seurat(seurat_rds)
        }
        
        if (params.notebook) {
            notebook = file("${projectDir}/notebooks/qc_plots.ipynb")
            run_notebook(notebook, cluster_seurat.out.rds)
        }
}