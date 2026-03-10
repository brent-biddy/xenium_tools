
include { SEURAT as SEURAT_OBJ } from "${projectDir}/modules/seurat.nf"
include { execute_notebook } from "${projectDir}/modules/notebook_execution.nf"


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
    def downsample_flag = params.downsample ? "--downsample" : ""
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

        if (params.run_create_bpcells) {
            create_bpcells_seurat_object(sample_info)
            bpcells_rds_ch = create_bpcells_seurat_object.out.full_rds
        } else {
            bpcells_rds_ch = Channel.empty()
        }

        if (params.run_qc_plots) {
            notebook_file = file("${projectDir}/notebooks/xenium_qc_plots.ipynb")
            execute_notebook(seurat_rds_ch, notebook_file, "xenium_qc_report")
        }

        if (params.run_create_seurat && params.run_create_bpcells) {
            joined = sample_info
                .join(seurat_rds_ch)
                .join(bpcells_rds_ch)
                .flatMap{ tuple(id: it[0], xenium_dir: it[1], seurat_rds: it[2], bpcells_rds: it[3]) }
            SEURAT_OBJ(joined)
        }
}