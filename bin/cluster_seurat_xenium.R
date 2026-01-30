#!/usr/bin/env Rscript

library(Seurat, quietly=TRUE)
library(data.table, quietly=TRUE)
library(dplyr, quietly=TRUE)
library(future, quietly=TRUE)

arguments <- function(){
    
    library(argparse)

    parser <- ArgumentParser(description='Create a Seurat object from Xenium data')
    parser$add_argument('--seurat_object', required=TRUE, help='Path to the Xenium data directory')
    parser$add_argument('--executor', required = FALSE, default = 'local', help = "Executor used when running with nextflow.")
    args <- parser$parse_args()
    return(args)    
}

cluster_seurat <- function(seurat_obj){
    print("Normalizing Data")
    seurat_obj <- NormalizeData(object =seurat_obj)
    print("Finding Variable Features")
    seurat_obj <- FindVariableFeatures(object =seurat_obj)
    print("Scaling Data")
    seurat_obj <- ScaleData(object =seurat_obj)
    print("Running PCA")
    seurat_obj <- RunPCA(object =seurat_obj)
    print("Finding Neighbors")
    seurat_obj <- FindNeighbors(object =seurat_obj, dims = 1:30)
    print("Clustering Data")
    seurat_obj <- FindClusters(object =seurat_obj)
    print("Running UMAP")
    seurat_obj <- RunUMAP(object =seurat_obj, dims = 1:30)
    print("UMAP Finished")
    return(seurat_obj)
}


main <- function(){

    args <- arguments()

    if(args$executor == 'local'){
        options(future.globals.maxSize = 2e9)
        plan(multisession, workers = availableCores(omit = 1))
        print(plan())
    } else if (args$executor == 'slurm') {
        options(future.globals.maxSize = 2e9)
        plan(multicore, workers = availableCores(omit = 1, methods = 'Slurm'))
        print(plan())
    }
    
    print("Loading Seurat Object")
    seurat <- readRDS(args$seurat_object)
    
    print("Beginning Clustering Workflow")
    seurat <- cluster_seurat(seurat)

    cluster_df <- seurat[[c("cell_id", "seurat_clusters")]]

    print("Writing Cluster info to CSV")
    write.csv(cluster_df, "seurat_clusters.csv", row.names = FALSE, quote = FALSE)

    print("Saving Seurat Object")
    saveRDS(seurat, "seurat_clusters.RDS")
}

main()
