#!/usr/bin/env Rscript

library(Seurat, quietly = TRUE)
library(data.table, quietly = TRUE)
library(argparse, quietly = TRUE)


arguments <- function(){

    parser <- ArgumentParser(description='Spatially-aware downsampling of a Xenium Seurat object')
    parser$add_argument('--rds_file', required=TRUE, help='Path to the input Seurat RDS file')
    parser$add_argument('--bin_size', required=FALSE, type='double', default=50,
                        help='Side length of spatial bins in microns (default: 50)')
    parser$add_argument('--fraction', required=FALSE, type='double', default=0.1,
                        help='Fraction of cells to sample per bin (default: 0.1)')

    args <- parser$parse_args()

    return(args)
}


spatial_downsample <- function(seurat_obj, bin_size = 100, fraction = 0.1){

    print(paste("Starting spatial downsampling with bin_size =", bin_size, "and fraction =", fraction))

    meta <- seurat_obj[[]]

    if(!all(c("x_centroid", "y_centroid") %in% colnames(meta))){
        stop("Seurat object metadata must contain 'x_centroid' and 'y_centroid' columns.")
    }

    meta$bin_x <- floor(meta$x_centroid / bin_size)
    meta$bin_y <- floor(meta$y_centroid / bin_size)
    meta$bin_id <- paste(meta$bin_x, meta$bin_y, sep = "_")
    meta$cell_id_rowname <- rownames(meta)

    bins <- split(meta$cell_id_rowname, meta$bin_id)

    sampled_cells <- unlist(lapply(bins, function(cells_in_bin){
        n_target <- max(1, floor(length(cells_in_bin) * fraction))
        if(length(cells_in_bin) <= n_target){
            return(cells_in_bin)
        }
        sample(cells_in_bin, n_target)
    }), use.names = FALSE)

    print(paste("Original cell count:", ncol(seurat_obj)))
    print(paste("Downsampled cell count:", length(sampled_cells)))
    print(paste("Total bins:", length(bins)))

    downsampled_obj <- subset(seurat_obj, cells = sampled_cells)

    return(downsampled_obj)
}


main <- function(){

    args <- arguments()

    print(paste("Loading Seurat object from:", args$rds_file))
    seurat_obj <- readRDS(args$rds_file)

    downsampled_obj <- spatial_downsample(seurat_obj, bin_size = args$bin_size, fraction = args$fraction)

    print("Saving downsampled Seurat object")
    saveRDS(downsampled_obj, file = "seurat_object_downsampled.RDS")
}


main()
