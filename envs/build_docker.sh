#!/usr/bin/env bash

mkdir images
docker build ./ --file ./Dockerfile.seurat --tag xenium_tools_seurat:latest