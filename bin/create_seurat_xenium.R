#!/usr/bin/env Rscript

library(Seurat, quietly=TRUE)
library(data.table, quietly=TRUE)
library(dplyr, quietly=TRUE)

arguments <- function(){
    
    library(argparse)

    parser <- ArgumentParser(description='Create a Seurat object from Xenium data')
    parser$add_argument('--data_dir', required=TRUE, help='Path to the Xenium data directory')
    parser$add_argument('--sample_name', required=FALSE, help='Sample name for the Seurat object', default='Xenium_Sample')
    parser$add_argument('--downsample', required = FALSE, action='store_true', help='Whether to downsample the data')
    args <- parser$parse_args()
    return(args)    
}

create_seurat_xenium <- function(data_dir, sample_name = "Xenium_Sample"){

    # output_file <- file.path(data_dir, paste0(sample_name, "_seurat.RDS"))
    output_file <- paste0(sample_name, "_seurat.RDS")
    cluster_file <- file.path(data_dir, "analysis", "clustering", "gene_expression_graphclust", "clusters.csv")
    umap_file <- file.path(data_dir, "analysis", "umap", "gene_expression_2_components", "projection.csv")
    cell_meta_file <- file.path(data_dir, "cells.csv.gz")

    print(paste("Creating Seurat object from Xenium data in: ", data_dir))
    print(paste("Output will be saved to: ", output_file))
    print(paste("Sample name set to: ", sample_name))

    # Load Xenium data
    seurat_obj <- LoadXenium(data.dir = data_dir, fov = "fov", segmentations = "cell")
    seurat_obj[["cell_id"]] <- rownames(seurat_obj[[]])
    seurat_obj$orig.ident <- sample_name
   
    if(file.exists(cluster_file)){
        print("Clustering data found. Adding Clustering Data to Seurat object.")
        print(paste0("Cluster data from the following file will be used: ", cluster_file))
        cluster_dat <- fread(cluster_file, data.table=FALSE, col.names=c("cell_id", "xenium_clusters"), colClasses=c("character", "integer"))
        seurat_obj[[]] <- left_join(seurat_obj[[]], cluster_dat, by="cell_id")
    }

    if(file.exists(cell_meta_file)){
        print("Cell Meta Data found. Adding meta data to Seurat object.")
        print(paste0("Meta data from the following file will be used: ", cell_meta_file))
        cell_meta_dat <- fread(cell_meta_file, data.table=FALSE, colClasses = c("cell_id"="character"))
        seurat_obj[[]] <- left_join(seurat_obj[[]][, !(colnames(seurat_obj[[]]) %in% c("segmentation_method"))], cell_meta_dat, by="cell_id")
    }

    if(file.exists(umap_file)){
        print("UMAP Projection found. Adding UMAP Projection to Seurat Object.")
        print(paste0("Projection from the following file will be used: ", umap_file))
        umap_dat <- fread(umap_file, col.names = c("Barcode", "UMAP_1", "UMAP_2"), data.table=FALSE)
        row.names(umap_dat) <- umap_dat[["Barcode"]]
        seurat_obj@reductions$xenium_umap <- SeuratObject::CreateDimReducObject(embeddings = as.matrix(umap_dat[c("UMAP_1", "UMAP_2")]), key = "UMAP_", assay = "Xenium")
    }
    
    if( file.exists(output_file) ){
        # output_file <- file.path(data_dir, paste0(sample_name, "_seurat_", as.integer(Sys.time()), ".RDS"))
        output_file <- paste0(sample_name, "_seurat_", as.integer(Sys.time()), ".RDS")
    }

    # Save Seurat object
    Idents(seurat_obj) <- "xenium_clusters"
    saveRDS(seurat_obj, file = "seurat_object.RDS")
    return(seurat_obj)
}

downsample_seurat_xenium <- function(seurat_obj, data_dir, sample_name="Xenium_Sample"){
    
    output_file <- file.path(data_dir, paste0(sample_name, "_seurat_downsampled.RDS"))

    print("Downsampling Seurat object")
    print(paste0("Downsampled Data will be saved to: ", output_file))
    
    seurat_obj <- subset(seurat_obj, downsample = 100)

    if( file.exists(output_file) ){
        # output_file <- file.path(data_dir, paste0(sample_name, "_seurat_downsampled_", as.integer(Sys.time()), ".RDS"))
        output_file <- paste0(sample_name, "_seurat_downsampled_", as.integer(Sys.time()), ".RDS")
    }

    saveRDS(seurat_obj, file = "seurat_object_downsampled.RDS")
}

main <- function(){

    args <- arguments()
    
    seurat <- create_seurat_xenium(data_dir = args$data_dir, sample_name = args$sample_name)

    if(args$downsample){
        seurat <- downsample_seurat_xenium(seurat_obj = seurat, data_dir = args$data_dir, sample_name = args$sample_name)
    }
}

main()
