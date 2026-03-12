process cluster_seurat {

    tag "${sample_name}"
    label 'seurat'

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

process sketch_cluster_seurat {

    tag "${sample_name}"
    label 'seurat'

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

workflow SEURAT {
    take:
        sample_info

    main:

        seurat_rds = sample_info.map{ tuple(it.id, it.seurat_rds) } // [sample_id, seurat_rds]
        bpcells_rds = sample_info.map{ tuple(it.id, it.bpcells_rds) } // [sample_id, bpcells_rds]
        xenium_dir = sample_info.map{ tuple(it.id, it.xenium_dir) } // [sample_id, xenium_dir]

        cluster_rds_ch = Channel.empty()

        if(params.cluster){
            sketch_cluster_seurat(seurat_rds)
        }

        if(params.cluster_full){
            cluster_seurat(seurat_rds)
            cluster_rds_ch = cluster_seurat.out.rds
        }

    emit:
        seurat_rds = seurat_rds
        cluster_rds = cluster_rds_ch
}