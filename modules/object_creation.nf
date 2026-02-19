
include { SEURAT as SEURAT_OBJ } from "${projectDir}/modules/seurat.nf"


process create_seurat_object {
    
    tag "${sample_name}"

    time = { 15.m * (1 + task.attempt)}

    input:
    tuple val(sample_name), path(xenium_output_path)

    output:
    tuple val(sample_name), path("seurat_object.RDS"), emit: "full_rds"
    tuple val(sample_name), path("seurat_object_downsampled.RDS"), emit: "small_rds", optional: true
    
    publishDir "${params.output_path}/results/${sample_name}", pattern: "seurat_object.RDS", saveAs: { "${sample_name}_seurat.RDS" }, mode: 'copy'
    publishDir "${params.output_path}/results/${sample_name}", pattern: "seurat_object_downsampled.RDS", saveAs: { "${sample_name}_seurat_downsampled.RDS" }, mode: 'copy'

    script:
    def downsample_flag = params.downsample == "true" ? "--downsample" : ""
    """
    create_seurat_xenium.R --data_dir ${xenium_output_path} --sample_name ${sample_name} ${downsample_flag}
    """

    stub:
    """
    touch seurat_object.RDS
    if [[ "${params.downsample}" == "true" ]]; then
        touch seurat_object_downsampled.RDS
    fi
    """
}

process create_bpcells_seurat_object {
    
    tag "${sample_name}"

    time = { 15.m * (1 + task.attempt)}

    input:
    tuple val(sample_name), path(xenium_output_path)

    output:
    tuple val(sample_name), path("bpcells_seurat_object.RDS"), emit: "full_rds"
    tuple val(sample_name), path("bpcells_seurat_object_downsampled.RDS"), emit: "small_rds", optional: true
    
    publishDir "${params.output_path}/results/${sample_name}/bpcells/", pattern: "bpcells_seurat_object.RDS", saveAs: { "${sample_name}_bpcells_seurat.RDS" }, mode: 'copy'
    publishDir "${params.output_path}/results/${sample_name}/bpcells/", pattern: "bpcells_seurat_object_downsampled.RDS", saveAs: { "${sample_name}_bpcells_seurat_downsampled.RDS" }, mode: 'copy'

    script:
    def downsample_flag = params.downsample == "true" ? "--downsample" : ""
    """
    create_bpcells_seurat_xenium.R --data_dir ${xenium_output_path} --sample_name ${sample_name} ${downsample_flag}
    """

    stub:
    """
    touch bpcells_seurat_object.RDS
    if [[ "${params.downsample}" == "true" ]]; then
        touch bpcells_seurat_object_downsampled.RDS
    fi
    """
}

process xenium_qc_plots {

    tag "${sample_name}"

    time = { 15.m * (1 + task.attempt)}

    input:
    path (notebook_path)
    tuple val(sample_name), path(seurat_rds)

    output:
    tuple val(sample_name), path("jupyter_notebook.html"), emit: html

    publishDir "${params.output_path}/results/${sample_name}", pattern: "jupyter_notebook.html", saveAs: { "${sample_name}_xenium_qc_report.html" }, mode: 'copy'

    script:
    """
    jupyter nbconvert --execute --allow-errors --output jupyter_notebook --to html ${notebook_path}
    """
    stub:
    """
    touch jupyter_notebook.html
    """
}

process bp_cells_clustering {

    tag "${sample_name}"

    time = { 15.m * (1 + task.attempt)}

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


workflow OBJECT_CREATION {
    take:
        sample_info
    main:
        
        create_seurat_object(sample_info)

        create_bpcells_seurat_object(sample_info)

        notebook_file = file("${projectDir}/notebooks/xenium_qc_plots.ipynb")
        xenium_qc_plots(notebook_file, create_seurat_object.out.full_rds)

        if(params.bp_clustering){
            bp_clustering_nb = file("${projectDir}/notebooks/bp_cells_clustering.ipynb")
            bp_cells_clustering(bp_clustering_nb, sample_info)
        }
        
        SEURAT_OBJ(create_seurat_object.out.full_rds)
}