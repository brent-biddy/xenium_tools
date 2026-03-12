
include { SEURAT as SEURAT_OBJ } from "${projectDir}/modules/seurat.nf"


process create_seurat_object {
    
    tag "${sample_name}"

    time = { 15.m * (1 + task.attempt)}

    input:
    tuple val(sample_name), path(xenium_output_path)

    output:
    tuple val(sample_name), path("seurat_object.RDS"), emit: "full_rds"

    publishDir "${params.output_path}/results/${sample_name}", pattern: "seurat_object.RDS", saveAs: { "${sample_name}_seurat.RDS" }, mode: 'copy'

    script:
    """
    create_seurat_xenium.R --data_dir ${xenium_output_path} --sample_name ${sample_name}
    """

    stub:
    """
    touch seurat_object.RDS
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
    def downsample_flag = params.downsample ? "--downsample" : ""
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


process downsample_seurat_object {

    tag "${sample_name}"

    time = { 15.m * (1 + task.attempt)}

    input:
    tuple val(sample_name), path(seurat_rds)

    output:
    tuple val(sample_name), path("seurat_object_downsampled.RDS"), emit: "downsampled_rds"

    publishDir "${params.output_path}/results/${sample_name}", pattern: "seurat_object_downsampled.RDS", saveAs: { "${sample_name}_seurat_downsampled.RDS" }, mode: 'copy'

    script:
    """
    downsample_seurat_xenium.R --rds_file ${seurat_rds} --bin_size ${params.downsample_bin_size} --fraction ${params.downsample_fraction}
    """

    stub:
    """
    touch seurat_object_downsampled.RDS
    """
}


process create_spatialdata_object {

    tag "${sample_name}"

    container 'babiddy755/xenium_tools_squidpy:latest'

    time = { 15.m * (1 + task.attempt)}

    input:
    tuple val(sample_name), path(xenium_output_path)

    output:
    tuple val(sample_name), path("data.zarr"), emit: "spatialdata_zarr"

    publishDir "${params.output_path}/results/${sample_name}", pattern: "data.zarr", saveAs: { "${sample_name}_spatialdata.zarr" }, mode: 'copy'

    script:
    """
    create_spatialdata_object.py ${xenium_output_path} --output_zarr data.zarr
    """

    stub:
    """
    mkdir -p data.zarr
    """
}


process downsample_spatialdata_object {

    tag "${sample_name}"

    container 'babiddy755/xenium_tools_squidpy:latest'

    time = { 15.m * (1 + task.attempt)}

    input:
    tuple val(sample_name), path(spatialdata_zarr)

    output:
    tuple val(sample_name), path("data_downsampled.zarr"), emit: "downsampled_zarr"

    publishDir "${params.output_path}/results/${sample_name}", pattern: "data_downsampled.zarr", saveAs: { "${sample_name}_spatialdata_downsampled.zarr" }, mode: 'copy'

    script:
    """
    downsample_spatialdata_object.py --zarr_file ${spatialdata_zarr} --bin_size ${params.downsample_bin_size} --fraction ${params.downsample_fraction}
    """

    stub:
    """
    mkdir -p data_downsampled.zarr
    """
}


workflow OBJECT_CREATION {
    take:
        sample_info
    main:
        if (params.run_create_seurat) {
            create_seurat_object(sample_info)
            seurat_rds_ch = create_seurat_object.out.full_rds
        } else {
            seurat_rds_ch = Channel.empty()
        }

        if (params.run_create_seurat && params.run_downsample) {
            downsample_seurat_object(seurat_rds_ch)
        }

        if (params.run_create_bpcells) {
            create_bpcells_seurat_object(sample_info)
            bpcells_rds_ch = create_bpcells_seurat_object.out.full_rds
        } else {
            bpcells_rds_ch = Channel.empty()
        }

        if (params.run_create_spatialdata) {
            create_spatialdata_object(sample_info)
            spatialdata_ch = create_spatialdata_object.out.spatialdata_zarr
        } else {
            spatialdata_ch = Channel.empty()
        }

        if (params.run_create_spatialdata && params.run_downsample) {
            downsample_spatialdata_object(spatialdata_ch)
        }

        if (params.run_create_seurat && params.run_create_bpcells) {
            joined = sample_info
                .join(seurat_rds_ch)
                .join(bpcells_rds_ch)
                .flatMap{ tuple(id: it[0], xenium_dir: it[1], seurat_rds: it[2], bpcells_rds: it[3]) }
            SEURAT_OBJ(joined)
        }

    emit:
        seurat_rds = seurat_rds_ch
        bpcells_rds = bpcells_rds_ch
        spatialdata = spatialdata_ch
}