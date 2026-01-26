#!/usr/bin/env Rscript

library(Seurat, quietly=TRUE)
library(data.table, quietly=TRUE)
library(dplyr, quietly=TRUE)

arguments <- function(){
    
    library(argparse)

    parser <- ArgumentParser(description='Create a Seurat object from Xenium data')
    parser$add_argument('--seurat_object', required=TRUE, help='Path to the Xenium data directory')
    parser$add_argument('--sample_name', required=FALSE, help='Sample name for the Seurat object', default='Xenium_Sample')
    args <- parser$parse_args()
    return(args)    
}

cluster_seurat <- function(seurat_obj, sample_name = "Xenium_Sample"){
    seurat_obj <- NormalizeData(object =seurat_obj)
    seurat_obj <- FindVariableFeatures(object =seurat_obj)
    seurat_obj <- ScaleData(object =seurat_obj)
    seurat_obj <- RunPCA(object =seurat_obj)
    seurat_obj <- FindNeighbors(object =seurat_obj, dims = 1:30)
    seurat_obj <- FindClusters(object =seurat_obj)
    seurat_obj <- RunUMAP(object =seurat_obj, dims = 1:30)
    return(seurat_obj)
}


main <- function(){

    args <- arguments()
    
    seurat <- readRDS(args$seurat_object)

    seurat <- cluster_seurat(seurat)

    saveRDS(seurat, "test_cluster_seurat.RDS")
}

main()
