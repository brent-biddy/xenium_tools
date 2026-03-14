# Container Build Instructions

This directory contains the container definition files for the Xenium analysis environment.

| File | Purpose |
|---|---|
| `Dockerfile.seurat` | Builds the seurat Docker image (R-based, used for QC, clustering, and marker scoring) |
| `Dockerfile.squidpy` | Builds the squidpy Docker image (Python-based, used for spatial analysis notebooks) |
| `squidpy_env.yml` | Conda environment definition for the squidpy container |

Throughout these instructions, replace `<dockerhub-username>` and `<image-name>` with your own Docker Hub username and a name of your choosing. After building and pushing, update the relevant `container` entry in `nextflow.config` to match your chosen image path.

---

## Prerequisites

- [Docker](https://docs.docker.com/get-docker/) installed and running
- A Docker Hub account with a repository to push to

---

## Building and Pushing Manually

### Seurat image

```bash
cd envs/
docker build -f Dockerfile.seurat -t <dockerhub-username>/<image-name>:latest .
```

> **Note:** If you hit GitHub API rate limits when installing R packages from GitHub during the build, you can pass a personal access token as a build secret. The token is only available during the build and is never baked into the image:
> ```bash
> GITHUB_PAT=<your-token> docker build \
>     --secret id=gh_token,env=GITHUB_PAT \
>     -f Dockerfile.seurat \
>     -t <dockerhub-username>/<image-name>:latest .
> ```

### Squidpy image

```bash
cd envs/
docker build -f Dockerfile.squidpy -t <dockerhub-username>/<image-name>:latest .
```

### Push to Docker Hub

```bash
docker login
docker push <dockerhub-username>/<image-name>:latest
```

### (Optional) Tag and push a versioned release

```bash
docker tag <dockerhub-username>/<image-name>:latest <dockerhub-username>/<image-name>:1.0.0
docker push <dockerhub-username>/<image-name>:1.0.0
```

---

## Running the Workflow

Once your images are built and `nextflow.config` is updated, run the workflow with your chosen profile:

```bash
# Docker
nextflow run main.nf --samplesheet samples.csv -profile docker

# Singularity (pulls from Docker Hub automatically)
nextflow run main.nf --samplesheet samples.csv -profile singularity

# Apptainer (pulls from Docker Hub automatically)
nextflow run main.nf --samplesheet samples.csv -profile apptainer

# OSCER / SLURM cluster
nextflow run main.nf --samplesheet samples.csv -profile oscer
```

---

## Container Contents

### Seurat

**Base Image:** `bioconductor/bioconductor_docker:3.21`

**System Dependencies:** `libhdf5-dev`, `libgeos-dev`, `libproj-dev`, `libudunits2-dev`, `libgdal-dev`, `libglpk-dev`, `pandoc`

**Python Packages:** `jupyter`, `nbconvert`, `ipykernel`, `matplotlib`

**R Packages:** `Seurat`, `SeuratObject`, `Signac`, `patchwork`, `data.table`, `dplyr`, `argparse`, `R.utils`, `arrow` (v21.0.0), `IRkernel`, `pals`

**Jupyter Kernel:** `r_seurat`

### Squidpy

**Base Image:** `mambaorg/micromamba:latest`

**Environment:** defined in `squidpy_env.yml` (conda-forge + bioconda)

**Packages:** `quarto`, `squidpy`, `spatialdata`, `spatialdata-io`, `scanpy`, `anndata`, `pandas`, `matplotlib`, `seaborn`, `jupyter`, `ipykernel`

**Jupyter Kernel:** `python_squidpy`

---

## Additional Resources

- [Docker Hub](https://hub.docker.com/)
- [Apptainer Documentation](https://apptainer.org/docs/)

---

## CI/CD: Manual Builds via GitHub Actions

Two workflows are available, each triggered manually from the GitHub Actions UI via the "Run workflow" button. Each builds the relevant Docker image, pushes it to Docker Hub tagged with both `latest` and the commit SHA, then automatically updates the corresponding `container` entry in `nextflow.config` and commits it back to the repo.

Layer caching is enabled via the GitHub Actions cache, so rebuilds that don't change early Dockerfile layers will be significantly faster.

| Workflow | File |
|---|---|
| Build and Push Seurat Docker Image | `.github/workflows/build-docker.yml` |
| Build and Push Squidpy Docker Image | `.github/workflows/build-docker-squidpy.yml` |

### Required GitHub Secrets and Variables

Before either action can run, configure the following in your repository under **Settings → Secrets and variables → Actions**:

| Type | Name | Used by |
|---|---|---|
| Secret | `DOCKERHUB_USERNAME` | Both |
| Secret | `DOCKERHUB_TOKEN` | Both |
| Secret | `GH_TOKEN` | Seurat only (avoids GitHub API rate limits for R package installs) |
| Variable | `DOCKERHUB_IMAGE` | Seurat (e.g. `xenium_tools_seurat`) |
| Variable | `DOCKERHUB_IMAGE_SQUIDPY` | Squidpy (e.g. `xenium_tools_squidpy`) |
