process cluster_seurat {
    
    tag "${sample_name}"

    time { 2.h * (1 + task.attempt)}
    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore'}

    input:
    tuple val(sample_name), path(seurat_obj)
    
    output:
    tuple val(sample_name), path("seurat_clusters.RDS"), emit: rds
    tuple val(sample_name), path("seurat_clusters.csv"), emit: csv
    
    publishDir "${params.output_path}/results/${sample_name}", pattern: "seurat_clusters.RDS", saveAs: { "${sample_name}_seurat_clusters.RDS" }, mode: 'copy'
    publishDir "${params.output_path}/results/${sample_name}", pattern: "seurat_clusters.csv", saveAs: { "${sample_name}_seurat_clusters.csv" }, mode: 'copy'


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

process seurat_cluster_plots {

    tag "${sample_name}"

    time = { 15.m * (1 + task.attempt) }

    input:
    path (notebook_path)
    tuple val(sample_name), path(seurat_rds)

    output:
    tuple val(sample_name), path("jupyter_notebook.html"), emit: html

    publishDir "${params.output_path}/results/${sample_name}", pattern: "jupyter_notebook.html", saveAs: { "${sample_name}_seurat_cluster_report.html" }, mode: 'copy'

    script:
    """
    jupyter nbconvert --execute --allow-errors --output jupyter_notebook --to html ${notebook_path}
    """
    stub:
    """
    touch jupyter_notebook.html
    """
}

process sketch_cluster_seurat {
    
    tag "${sample_name}"

    time = { 15.m * (1 + task.attempt)}

    input:
    tuple val(sample_name), path(seurat_obj)
    
    output:
    tuple val(sample_name), path("test_sketch_and_cluster_seurat.RDS"), emit: rds
    tuple val(sample_name), path("test_sketch_and_cluster_seurat_clusters.csv"), emit: csv
    
    publishDir "${params.output_path}/results/${sample_name}", pattern: "test_sketch_and_cluster_seurat.RDS", saveAs: { "${sample_name}_seurat_sketch_clusters.RDS" }, mode: 'copy'
    publishDir "${params.output_path}/results/${sample_name}", pattern: "test_sketch_and_cluster_seurat_clusters.csv", saveAs: { "${sample_name}_seurat_sketch_clusters.csv" }, mode: 'copy'


    script:
    """
    sketch_cluster_xenium.R --seurat_object ${seurat_obj}
    """
    stub:
    """
    touch test_sketch_and_cluster_seurat.RDS
    touch test_sketch_and_cluster_seurat_clusters.csv
    """
}

process seurat_score_markers {

    tag "${sample_name}"

    time = { 15.m * (1 + task.attempt)}

    input:
    path (notebook_path)
    path(yaml_path)
    tuple val(sample_name), path(seurat_rds)

    output:
    tuple val(sample_name), path("jupyter_notebook.html"), emit: html

    publishDir "${params.output_path}/results/${sample_name}", pattern: "jupyter_notebook.html", saveAs: { "${sample_name}_celltype_marker_gene_report.html" }, mode: 'copy'

    script:
    """
    jupyter nbconvert --execute --allow-errors --output jupyter_notebook --to html ${notebook_path}
    """
    stub:
    """
    touch jupyter_notebook.html
    """
}

process seurat_subcluster_notebook{

    tag "${sample_name}"
    
    time = { 15.m * (1 + task.attempt)}
    
    input:
    path (notebook_path)
    path(yaml_path)
    path(annotation_csv)
    tuple val(sample_name), path(seurat_rds)

    output:
    tuple val(sample_name), path("jupyter_notebook.html"), emit: html

    publishDir "${params.output_path}/results/notebooks/${sample_name}", pattern: "jupyter_notebook.html", mode: 'copy'

    script:
    """
    Rscript -e 'rmarkdown::render(input = "${notebook_path}",\
                                    params = list(marker_path = "${yaml_path}",\
                                    annotation_path = "${annotation_csv}",\
                                    rds_path = "${seurat_rds}")\
                                )'
    """
    stub:
    """
    touch jupyter_notebook.html
    """

}


workflow SEURAT {
    take:
        seurat_rds

    main:
        if(params.cluster){
            sketch_cluster_seurat(seurat_rds)
        }
        
        if(params.cluster_full){
            cluster_seurat(seurat_rds)
            cluster_notebook = file("${projectDir}/notebooks/seurat_cluster_plots.ipynb")
            seurat_cluster_plots(cluster_notebook, cluster_seurat.out.rds)
        }

        if(params.score_markers){
            marker_notebook = file("${projectDir}/notebooks/marker_scores.ipynb")
            marker_yaml = file("${projectDir}/refs/ovary_markers.yaml")
            seurat_score_markers(marker_notebook, marker_yaml, seurat_rds)
        }

        if(params.subcluster_nb){
            subcluster_notebook = file("${projectDir}/notebooks/TEST_roi_3_marker_scores_subcluster_oocyte_seurat_clusters.Rmd")
            annotation_file = file("${projectDir}/refs/follicle_annotations.csv")
            marker_yaml = file("${projectDir}/refs/ovary_markers.yaml")
            seurat_subcluster_notebook(subcluster_notebook, marker_yaml, annotation_file, seurat_rds)
        }
}