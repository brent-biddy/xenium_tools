// main.nf

process build_seurat_sif {
    
    tag "build_seurat_sif"

    input:
    path seurat_def

    output:
    path "seurat.sif"

    publishDir "${projectDir}/envs/images", mode: 'copy'

    script:
    """
    singularity build --fakeroot seurat.sif ${seurat_def}
    """
    stub:
    """
    touch seurat.sif
    """
}

process create_seurat_object {
    
    tag "${sample_name}"

    container {workflow.containerEngine == 'singularity' ? // If using Singularity
                "${projectDir}/envs/images/seurat.sif" : //Then container to use is this line.
                "babiddy755/xenium_tools_seurat:latest"}// Else, container to use is this line.

    input:
    tuple val(sample_name), path(xenium_output_path), val(downsample)
    
    output:
    tuple val(sample_name), path("seurat_object.RDS")
    tuple val(sample_name), path("seurat_object_downsampled.RDS"), optional: true
    
    publishDir "results/${sample_name}", mode: 'copy'

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

    container {workflow.containerEngine == 'singularity' ? // If using Singularity
                "${projectDir}/envs/images/seurat.sif" : //Then container to use is this line.
                "babiddy755/xenium_tools_seurat:latest"}// Else, container to use is this line.

    input:
    tuple val(sample_name), path(seurat_obj)
    
    output:
    tuple val(sample_name), path("test_seurat_cluster.RDS")
    
    publishDir "results/${sample_name}", mode: 'copy'

    script:
    """
    cluster_seurat_xenium.R --seurat_object ${seurat_obj} --sample_name ${sample_name}
    """
    stub:
    """
    touch test_cluster_seurat.RDS
    """
}

// Workflow block
workflow {

    if (params.build_seurat) {
        Channel.fromPath("${projectDir}/envs/seurat.def") |
        set{ seurat_def }
        build_seurat_sif(seurat_def)
    } else {

        Channel.fromPath(params.samplesheet) |
        splitCsv(header:true) |
        map{ row -> tuple(row.sample, file(row.path), params.downsample) } |
        set{sample_info}
        sample_info.dump(tag: "sample_info")

        create_seurat_object(sample_info)
    }
}