
include { SEURAT as SEURAT_OBJ } from "${projectDir}/modules/seurat.nf"


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

process xenium_qc_plots {

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


workflow OBJECT_CREATION {
    take:
        sample_info
    main:
        
        create_seurat_object(sample_info, params.downsample)

        notebook_file = file("${projectDir}/notebooks/xenium_qc_plots.ipynb")
        xenium_qc_plots(notebook_file, create_seurat_object.out.full_rds)
        
        SEURAT_OBJ(create_seurat_object.out.full_rds)
}