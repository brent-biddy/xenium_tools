// main.nf

process create_seurat_object {
    
    tag "${sample_name}"

    input:
    tuple val(sample_name), path(xenium_output_path), val(downsample)
    
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
    tuple val(sample_name), path("seurat_clusters.RDS")
    tuple val(sample_name), path("seurat_clusters.csv")
    
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

// Workflow block
workflow {

    Channel.fromPath(params.samplesheet) |
    splitCsv(header:true) |
    map{ row -> tuple(row.sample, file(row.path), params.downsample) } |
    set{sample_info}

    sample_info.dump(tag: "sample_info")
    
    create_seurat_object(sample_info)

    if (params.cluster) {
        cluster_seurat(create_seurat_object.out.full_rds)
    }
}