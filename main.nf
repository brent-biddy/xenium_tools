// main.nf

include { OBJECT_CREATION } from "${projectDir}/modules/object_creation.nf"
include { SEURAT as SEURAT_RDS } from "${projectDir}/modules/seurat.nf"


// Workflow block
workflow {

    if ( ! params.samplesheet ){
        error "ERROR: --samplesheet parameter is required"
    }

    if ( !file(params.samplesheet).exists()){
        error "ERROR: Samplesheet file not found: ${params.samplesheet}"
    }

    if (params.run_qc_plots && !params.run_create_seurat) {
        error "ERROR: run_qc_plots requires run_create_seurat to be true"
    }

    if (params.run_cluster_plots && !params.cluster_full) {
        log.warn "WARNING: run_cluster_plots=true has no effect when cluster_full=false"
    }

    if (params.score_markers && !file(params.marker_yaml).exists()) {
        error "ERROR: score_markers is enabled but marker_yaml file not found: ${params.marker_yaml}"
    }

    Channel.fromPath(params.samplesheet)
        .splitCsv(header:true)
        .map{ row -> 
            if (!row.sample || !row.path) { 
                error "ERROR: Samplesheet must have 'sample' and 'path' columns"
            }
            tuple(row.sample, file(row.path)) 
        }
        .branch{
            seurat: it[1].getExtension() == "RDS"
            anndata: it[1].getExtension() == "h5ad"
            dir: it[1].isDirectory()
        }
    .set{sample_info}

    sample_info.seurat.dump(tag: "seurat_input")
    sample_info.anndata.dump(tag: "anndata_input")
    sample_info.dir.dump(tag: "dir_input")

    if (params.run_object_creation) {
        OBJECT_CREATION(sample_info.dir)
    }

    if (params.run_seurat_rds) {
        SEURAT_RDS(sample_info.seurat)
    }

}

workflow.onComplete {
    println """
    Pipeline execution summary
    ---------------------------
    Completed at : ${workflow.complete}
    Duration     : ${workflow.duration}
    Success      : ${workflow.success}
    Work Dir     : ${workflow.workDir}
    Results Dir  : ${params.output_path}
    Exit status  : ${workflow.exitStatus}
    Error report : ${workflow.errorReport ?: '-'}
    """.stripIndent()
}

workflow.onError {
    println "Oops... Pipeline execution stopped with the following message: ${workflow.errorMessage}"
}

