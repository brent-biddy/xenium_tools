#!/usr/bin/env Rscript

arguments <- function(){
    
    library(argparse, quietly = TRUE)

    parser <- ArgumentParser(description='Create a Seurat object from Xenium data')
    parser$add_argument('--seurat_object', required=TRUE, help='Path to the Xenium data directory')
    parser$add_argument('--sample_name', required=FALSE, help='Sample name for the Seurat object', default='Xenium_Sample')
    args <- parser$parse_args()
    return(args)    
}

sketch_and_cluster_seurat <- function(seurat_obj, sample_name = "Xenium_Sample", sketch_method = "Uniform"){
    Idents(seurat_obj) <- "xenium_clusters"
    seurat_obj <- NormalizeData(object =seurat_obj)
    seurat_obj <- FindVariableFeatures(object =seurat_obj)
    seurat_obj <- SketchData(object = seurat_obj, ncells = 50000, method = "Uniform", sketched.assay = "sketch")
    DefaultAssay(seurat_obj) <- "sketch"
    seurat_obj <- FindVariableFeatures(object =seurat_obj)
    seurat_obj <- ScaleData(object =seurat_obj)
    seurat_obj <- RunPCA(object =seurat_obj, reduction.name = "pca.sketch", reduction.key = "pca_")
    seurat_obj <- FindNeighbors(object =seurat_obj, dims = 1:30, reduction = "pca.sketch")
    seurat_obj <- FindClusters(object =seurat_obj)
    seurat_obj <- RunUMAP(object =seurat_obj, dims = 1:30, reduction = "pca.sketch", return.model = TRUE)
    seurat_obj <- ProjectData(object = seurat_obj, assa = "Xenium", full.reduction = "pca.full", sketched.assay = "sketch", sketched.reduction = "pca.sketch", umap.model = "umap", dims = 1:30, refdata = list(projected_clusters = "seurat_clusters"))
    DefaultAssay(seurat_obj) <- "Xenium"
    return(seurat_obj)
}


main <- function(){

    args <- arguments()
    
    #Load Libraries after checking arguments.
    library(Seurat, quietly=TRUE)
    library(data.table, quietly=TRUE)
    library(dplyr, quietly=TRUE)
    options(future.globals.maxSize = 2e9)
    
    seurat <- readRDS(args$seurat_object)

    seurat <- sketch_and_cluster_seurat(seurat)

    cluster_df <- seurat[[c("cell_id", "projected_clusters")]]

    print("Writing Cluster infor to CSV")
    write.csv(cluster_df, "test_sketch_and_cluster_seurat_clusters.csv", row.names = FALSE, quote = FALSE)

    saveRDS(seurat, "test_sketch_and_cluster_seurat.RDS")
}

main()
