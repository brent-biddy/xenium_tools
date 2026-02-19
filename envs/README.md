# Container Build Instructions

This directory contains the container definition files for the Xenium analysis environment.

| File | Purpose |
|---|---|
| `Dockerfile.seurat` | Builds a Docker image for push to Docker Hub (used by all Nextflow profiles) |
| `seurat.def` | Builds a Singularity `.sif` image for push to Sylabs Cloud (optional alternative) |

Throughout these instructions, replace `<dockerhub-username>` and `<image-name>` with your own Docker Hub username and a name of your choosing, and replace `<sylabs-username>` and `<collection>` with your own Sylabs credentials. After building and pushing, update `process.container` in `nextflow.config` to match your chosen image path.

---

## Prerequisites

- [Docker](https://docs.docker.com/get-docker/) installed and running (for Docker builds)
- [Singularity](https://sylabs.io/docs/) or [Apptainer](https://apptainer.org/docs/) v3.0+ (for `.sif` builds)
- A Docker Hub account with a repository to push to
- (Optional) A Sylabs Cloud account if pushing `.sif` images

---

## Option 1: Docker Hub (Recommended)

This is the primary method. All Nextflow profiles (`docker`, `singularity`, `apptainer`, `oscer`) pull from Docker Hub. Singularity and Apptainer can pull Docker Hub images natively, so a single Docker build covers all execution environments.

### 1. Build the Docker image

```bash
cd envs/
docker build -f Dockerfile.seurat -t <dockerhub-username>/<image-name>:latest .
```

### 2. Log in to Docker Hub

```bash
docker login
```

### 3. Push the image to Docker Hub

```bash
docker push <dockerhub-username>/<image-name>:latest
```

### 4. Update `nextflow.config`

Set `process.container` to your image path in each relevant profile:

```groovy
process.container = '<dockerhub-username>/<image-name>:latest'
```

### 5. (Optional) Tag and push a versioned release

```bash
docker tag <dockerhub-username>/<image-name>:latest <dockerhub-username>/<image-name>:1.0.0
docker push <dockerhub-username>/<image-name>:1.0.0
```

---

## Option 2: Sylabs Cloud (Singularity/Apptainer `.sif`)

Use this if you prefer to distribute a `.sif` image via Sylabs Cloud rather than Docker Hub.

### 1. Log in to Sylabs Cloud

```bash
singularity remote login
# Paste your access token from https://cloud.sylabs.io/auth/tokens when prompted
```

### 2. Build the `.sif` image

```bash
cd envs/
mkdir -p images
# If you have root access, --fakeroot can be replaced with sudo
singularity build --fakeroot images/seurat.sif seurat.def
```

### 3. (Optional) Sign the image

```bash
singularity sign images/seurat.sif
```

### 4. Push to Sylabs Cloud

```bash
singularity push images/seurat.sif library://<sylabs-username>/<collection>/seurat:latest
```

### 5. Update `nextflow.config`

Set `process.container` to your Sylabs library path in each relevant profile:

```groovy
process.container = 'library://<sylabs-username>/<collection>/seurat:latest'
```

---

## Option 3: Local Singularity `.sif` File (HPC without Registry Access)

This option is intended for HPC environments where Docker is not available and outbound internet access is restricted or pulling from a registry is not permitted. Build the `.sif` image on a machine with internet access using the steps in Option 2, transfer the file to your HPC, and reference it by path.

After building and transferring `images/seurat.sif`, edit the singularity profile in `nextflow.config` to point to the local path:

```groovy
singularity {
    singularity.enabled = true
    process.container = '/full/path/to/envs/images/seurat.sif'
    process.executor = 'local'
}
```

> **Tip:** Use an absolute path for the `.sif` file — relative paths can cause issues depending on where Nextflow is invoked from.

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
- [Singularity Documentation](https://sylabs.io/docs/)
- [Sylabs Cloud](https://cloud.sylabs.io/)
- [Apptainer Documentation](https://apptainer.org/docs/)

---

---

## CI/CD: Manual Builds via GitHub Actions

The workflow at `.github/workflows/build-docker.yml` can be triggered manually from the GitHub Actions UI via the "Run workflow" button. It builds, verifies, and pushes the Docker image to Docker Hub. Package verification runs inside the built container before any push occurs, so a broken build will never reach Docker Hub.

### Required GitHub Secrets and Variables

Before the action can run, configure the following in your repository under **Settings → Secrets and variables → Actions**:

| Type | Name | Value |
|---|---|---|
| Secret | `DOCKERHUB_USERNAME` | Your Docker Hub username |
| Secret | `DOCKERHUB_TOKEN` | A Docker Hub [access token](https://app.docker.com/settings/personal-access-tokens) (not your password) |
| Variable | `DOCKERHUB_IMAGE` | The image name to push to (e.g. `xenium_tools_seurat`) |

Once set, trigger the workflow from **Actions → Build and Push Docker Image → Run workflow**. The image is tagged with both `latest` and the commit SHA, and `nextflow.config` is automatically updated with the new image reference and committed back to the repo.
