#!/usr/bin/env Rscript

library(Seurat, quietly = TRUE)
library(data.table, quietly = TRUE)
library(dplyr, quietly = TRUE)
library(argparse, quietly = TRUE)
library(SeuratObject, quietly = TRUE)


arguments <- function(){
    
    parser <- ArgumentParser(description='Create a Seurat object from Xenium data')
    parser$add_argument('--data_dir', required=TRUE, help='Path to the Xenium data directory')
    parser$add_argument('--sample_name', required=FALSE, help='Sample name for the Seurat object', default='Xenium_Sample')
    parser$add_argument('--downsample', required = FALSE, action='store_true', help='Whether to downsample the data')

    args <- parser$parse_args()
    
    return(args)    
}

add_meta_and_clusters <- function(seurat_obj, data_dir){

    print("Adding Cell Meta Data and Cluster Information to Seurat Object")
    cell_meta_file <- file.path(data_dir, "cells.csv.gz")
    cluster_file <- file.path(data_dir, "analysis", "clustering", "gene_expression_graphclust", "clusters.csv")
    
    if(file.exists(cluster_file)){

        print("Clustering data found. Adding Clustering Data to Seurat object.")
        print(paste("Cluster data from the following file will be used: ", cluster_file))
        cluster_dat <- fread(cluster_file, data.table=FALSE, col.names=c("cell_id", "xenium_clusters"), colClasses=c("character", "integer"))
        seurat_obj[[]] <- left_join(seurat_obj[[]], cluster_dat, by="cell_id")
    }

    if(file.exists(cell_meta_file)){

        print("Cell Meta Data found. Adding meta data to Seurat object.")
        print(paste("Meta data from the following file will be used: ", cell_meta_file))
        cell_meta_dat <- fread(cell_meta_file, data.table=FALSE, colClasses = c("cell_id"="character"))
        seurat_obj[[]] <- left_join(seurat_obj[[]][, !(colnames(seurat_obj[[]]) %in% c("segmentation_method"))], cell_meta_dat, by="cell_id")
    }

    Idents(seurat_obj) <- "xenium_clusters"

    return(seurat_obj)

}

add_reductions <- function(seurat_obj, data_dir){
    
    print("Adding Reduction to Seurat Object")
    pca_dir <- file.path(data_dir, "analysis", "pca", "gene_expression_10_components")
    umap_file <- file.path(data_dir, "analysis", "umap", "gene_expression_2_components", "projection.csv")
    
    if(dir.exists(pca_dir)){

        print(paste("Adding PCA to the Seurat Object from the following directory: ", pca_dir))        
        embedding_path <- file.path(pca_dir, "projection.csv")
        loadings_path <- file.path(pca_dir, "components.csv")

        embeddings <- fread(embedding_path, data.table = FALSE)
        colnames(embeddings) <- gsub(pattern = "-", replacement = "_", x = colnames(embeddings))
        row.names(embeddings) <- embeddings$Barcode
        embeddings <- as.matrix(embeddings[, 2:ncol(embeddings)])

        loadings <- fread(loadings_path, data.table = FALSE)
        loadings$PC <- paste0("PC_", loadings$PC)
        row.names(loadings) <- loadings$PC
        loadings <- t(as.matrix(loadings[, 2:ncol(loadings)]))

        pca <- CreateDimReducObject(embeddings = embeddings, loadings = loadings, assay = "Xenium", key = "PC_")

        seurat_obj@reductions$xenium_pca <- pca
    }
    
    if(file.exists(umap_file)){

        print(paste("Adding UMAP to Seurat Object from the following file: ", umap_file))
        umap_dat <- fread(umap_file, data.table=FALSE)
        colnames(umap_dat) <- gsub(pattern = "-", replacement = "_", x = colnames(umap_dat))
        row.names(umap_dat) <- umap_dat[["Barcode"]]
        umap_dat <- as.matrix(umap_dat[, 2:ncol(umap_dat)])

        umap <- CreateDimReducObject(embeddings = umap_dat, key = "UMAP_", assay = "Xenium")
        seurat_obj@reductions$xenium_umap <- umap
    }

    return(seurat_obj)

}

create_seurat_xenium <- function(data_dir, sample_name = "Xenium_Sample", downsample = TRUE){

    print(paste("Creating Seurat object from Xenium data in: ", data_dir))
    print(paste("Sample name set to: ", sample_name))

    seurat_obj <- LoadXenium(data.dir = data_dir, fov = "fov", segmentations = "cell")
    seurat_obj[["cell_id"]] <- rownames(seurat_obj[[]])
    seurat_obj$orig.ident <- sample_name

    seurat_obj <- add_meta_and_clusters(seurat_obj, data_dir)
    seurat_obj <- add_reductions(seurat_obj, data_dir)

    if(downsample){
        print("Saving Downsampled Seurat Object")
        saveRDS(subset(seurat_obj, downsample = 100), file = "seurat_object_downsampled.RDS")
    }

    print("Saving Full Seurat Object")
    saveRDS(seurat_obj, file = "seurat_object.RDS")
}

main <- function(){

    args <- arguments()
    
    seurat <- create_seurat_xenium(data_dir = args$data_dir, sample_name = args$sample_name)
}

main()
