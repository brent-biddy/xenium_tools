#!/usr/bin/env Rscript

library(Seurat, quietly=TRUE)
library(data.table, quietly=TRUE)
library(dplyr, quietly=TRUE)

arguments <- function(){
    
    library(argparse)

    parser <- ArgumentParser(description='Create a Seurat object from Xenium data')
    parser$add_argument('--seurat_object', required=TRUE, help='Path to the Xenium data directory')
    args <- parser$parse_args()
    return(args)    
}

cluster_seurat <- function(seurat_obj){
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

    cluster_df <- seurat[[c("cell_id", "seurat_clusters")]]

    write.csv(cluster_df, "seurat_clusters.csv", row.names = FALSE, quote = FALSE)

    saveRDS(seurat, "seurat_clusters.RDS")
}

main()
