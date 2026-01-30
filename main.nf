// main.nf

include { OBJECT_CREATION } from "${projectDir}/modules/object_creation.nf"
include { SEURAT as SEURAT_RDS } from "${projectDir}/modules/seurat.nf"


// Workflow block
workflow {

    Channel.fromPath(params.samplesheet)
        .splitCsv(header:true)
        .map{ row -> tuple(row.sample, file(row.path)) }
        .branch{
            seurat: it[1].getExtension() == "RDS"
            anndata: it[1].getExtension() == "h5ad"
            dir: it[1].isDirectory()
        }
        .set{sample_info}

    sample_info.seurat.dump(tag: "seurat_input")
    sample_info.anndata.dump(tag: "anndata_input")
    sample_info.dir.dump(tag: "dir_input")

    OBJECT_CREATION(sample_info.dir)
    SEURAT_RDS(sample_info.seurat)

}