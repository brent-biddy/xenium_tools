# Container Build Instructions

This directory contains the container definition files for the Xenium analysis environment.

| File | Purpose |
|---|---|
| `Dockerfile.seurat` | Builds a Docker image for push to Docker Hub (used by all Nextflow profiles) |

Throughout these instructions, replace `<dockerhub-username>` and `<image-name>` with your own Docker Hub username and a name of your choosing. After building and pushing, update `process.container` in `nextflow.config` to match your chosen image path.

---

## Prerequisites

- [Docker](https://docs.docker.com/get-docker/) installed and running
- A Docker Hub account with a repository to push to

---

## Building and Pushing Manually

### 1. Build the Docker image

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

### 2. Log in to Docker Hub

```bash
docker login
```

### 3. Push the image to Docker Hub

```bash
docker push <dockerhub-username>/<image-name>:latest
```

### 4. Update `nextflow.config`

Set `process.container` to your image path:

```groovy
process.container = '<dockerhub-username>/<image-name>:latest'
```

### 5. (Optional) Tag and push a versioned release

```bash
docker tag <dockerhub-username>/<image-name>:latest <dockerhub-username>/<image-name>:1.0.0
docker push <dockerhub-username>/<image-name>:1.0.0
```

---

## Running the Workflow

Once your image is built and `nextflow.config` is updated, run the workflow with your chosen profile:

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

**Base Image:** `bioconductor/bioconductor_docker:3.21`

**System Dependencies:** `libhdf5-dev`, `libgeos-dev`, `libproj-dev`, `libudunits2-dev`, `libgdal-dev`, `libglpk-dev`, `pandoc`

**Python Packages:** `jupyter`, `nbconvert`, `ipykernel`, `matplotlib`

**R Packages:** `Seurat`, `SeuratObject`, `Signac`, `patchwork`, `data.table`, `dplyr`, `argparse`, `R.utils`, `arrow` (v21.0.0), `IRkernel`, `pals`

---

## Additional Resources

- [Docker Hub](https://hub.docker.com/)
- [Apptainer Documentation](https://apptainer.org/docs/)

---

## CI/CD: Manual Builds via GitHub Actions

The workflow at `.github/workflows/build-docker.yml` is triggered manually from the GitHub Actions UI via the "Run workflow" button. It builds the Docker image and pushes it to Docker Hub, tagging it with both `latest` and the commit SHA. After a successful push, `nextflow.config` is automatically updated with the new SHA-tagged image reference and committed back to the repo.

Layer caching is enabled via the GitHub Actions cache, so rebuilds that don't change early Dockerfile layers (base image, system dependencies) will be significantly faster.

### Required GitHub Secrets and Variables

Before the action can run, configure the following in your repository under **Settings → Secrets and variables → Actions**:

| Type | Name | Value |
|---|---|---|
| Secret | `DOCKERHUB_USERNAME` | Your Docker Hub username |
| Secret | `DOCKERHUB_TOKEN` | A Docker Hub [access token](https://app.docker.com/settings/personal-access-tokens) (not your password) |
| Secret | `GH_TOKEN` | A GitHub [personal access token](https://github.com/settings/tokens) (optional — used to avoid GitHub API rate limits when installing R packages from GitHub during the build) |
| Variable | `DOCKERHUB_IMAGE` | The image name to push to (e.g. `xenium_tools_seurat`) |

Once set, trigger the workflow from **Actions → Build and Push Docker Image → Run workflow**.
