#!/usr/bin/env bash

mkdir images
singularity build --fakeroot ./images/seurat.sif ./seurat.def