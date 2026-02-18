# Container Build Instructions

This directory contains the Singularity/Apptainer definition file for the Xenium analysis environment.

## Prerequisites

- Singularity/Apptainer installed (v3.0+)
- `sudo` or `--fakeroot` capabilities for building
- (Optional) Sylabs Cloud account for pushing containers

## Quick Start: Build Locally

Build the container image:

```bash
cd envs/
mkdir -p images
singularity build --fakeroot images/seurat.sif seurat.def
```

## Use with Nextflow (Local)

After building locally, update `nextflow.config` to use your local container:

```groovy
profiles {
    singularity {
        docker.enabled = false
        singularity.enabled = true
        apptainer.enabled = false
        process.container = '/full/path/to/envs/images/seurat.sif'
        process.executor = 'local'
    }
}
```

Then run the workflow:
```bash
nextflow run main.nf --samplesheet samples.csv -profile singularity
```

## Push to Sylabs Cloud (Optional)

To share the container or use it across multiple systems, push to Sylabs Cloud.

**Prerequisites:**
- Sylabs Cloud account: https://cloud.sylabs.io/
- Access token generated from: https://cloud.sylabs.io/auth/tokens

### 1. Login to Sylabs Cloud

```bash
singularity remote login
# Paste your access token when prompted
```

### 2. Build and Sign the Container

```bash
# Build the container
mkdir -p images
singularity build images/seurat.sif seurat.def

# Sign the container (optional but recommended)
singularity sign images/seurat.sif
```

### 3. Push to Sylabs Cloud

```bash
# Push to your library
singularity push images/seurat.sif library://YOUR_USERNAME/xenium_tools/seurat:latest

# Example:
# singularity push images/seurat.sif library://babiddy/xenium_tools/seurat:latest
```

### 4. Update Nextflow Config

Update `nextflow.config` to use your Sylabs Cloud container:

```groovy
profiles {
    singularity {
        docker.enabled = false
        singularity.enabled = true
        apptainer.enabled = false
        process.container = 'library://YOUR_USERNAME/xenium_tools/seurat:latest'
        process.executor = 'local'
    }
}
```

## Pull and Use Remote Container

Anyone can now pull your container:

```bash
# Pull from Sylabs Cloud
singularity pull library://YOUR_USERNAME/xenium_tools/seurat:latest

# Or let Nextflow pull automatically when running the workflow
nextflow run main.nf --samplesheet samples.csv -profile singularity
```

## Container Contents

**Base Image:** bioconductor/bioconductor_docker:3.21

**System Dependencies:**
- libhdf5-dev, libgeos-dev, libproj-dev
- libudunits2-dev, libgdal-dev, libglpk-dev
- pandoc

**Python Packages:**
- jupyter, nbconvert, ipykernel, matplotlib

**R Packages:**
- Seurat, SeuratObject, Signac
- patchwork, data.table, dplyr
- argparse, R.utils
- arrow (v21.0.0)
- IRkernel

## Additional Resources

- [Singularity Documentation](https://sylabs.io/docs/)
- [Sylabs Cloud](https://cloud.sylabs.io/)
- [Apptainer Documentation](https://apptainer.org/docs/)
