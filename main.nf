// main.nf

process create_seurat_object {
    
    tag "${sample_name}"

    input:
    tuple val(sample_name), path(xenium_output_path)
    val(downsample)
    
    output:
    tuple val(sample_name), path("seurat_object.RDS"), emit: "full_rds"
    tuple val(sample_name), path("seurat_object_downsampled.RDS"), emit: "small_rds", optional: true
    
    publishDir "${params.output_path}/results/${sample_name}", pattern: "seurat_object.RDS", saveAs: { "${sample_name}_seurat.RDS" }, mode: 'copy'
    publishDir "${params.output_path}/results/${sample_name}", pattern: "seurat_object_downsampled.RDS", saveAs: { "${sample_name}_seurat_downsampled.RDS" }, mode: 'copy'

    script:
    def downsample = downsample ? "--downsample" : ""
    """
    create_seurat_xenium.R --data_dir ${xenium_output_path} --sample_name ${sample_name} ${downsample}
    """
    stub:
    """
    touch seurat_object.RDS
    if [[ $params.downsample ]]; then
        touch seurat_object_downsampled.RDS
    fi
    """
}

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
    cluster_seurat_xenium.R --seurat_object ${seurat_obj}
    """
    stub:
    """
    touch seurat_clusters.RDS
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

// Workflow block
workflow {

    Channel.fromPath(params.samplesheet)
        .splitCsv(header:true)
        .map{ row -> tuple(row.sample, file(row.path)) }
        .set{sample_info}

    sample_info.dump(tag: "sample_info")

    create_seurat_object(sample_info, params.downsample)

    if (params.cluster) {
        cluster_seurat(create_seurat_object.out.full_rds)
    }

    if (params.notebook) {
        notebook = file("${projectDir}/notebooks/test.ipynb")
        run_notebook(notebook, create_seurat_object.out.full_rds)
    }
}