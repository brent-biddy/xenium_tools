# xenium_tools

A Nextflow workflow for processing 10x Genomics Xenium spatial transcriptomics datasets. Accepts raw Xenium output directories or pre-built Seurat `.RDS` objects and produces clustered Seurat objects and QC/analysis reports.

## Requirements

- Nextflow ≥ 21.04.0
- Docker, Singularity, or Apptainer (see [Profiles](#profiles))

## Input

The workflow takes a CSV samplesheet with two columns:

| column | description |
|--------|-------------|
| `sample` | Sample name |
| `path` | Path to a Xenium output directory, a Seurat `.RDS` file, or an AnnData `.h5ad` file |

```csv
sample,path
sample1,/data/xenium/sample1
sample2,/data/seurat/sample2.RDS
```

## Usage

```bash
nextflow run main.nf \
  --samplesheet samplesheet.csv \
  --output_path ./results
```

To receive an email notification on completion or failure:

```bash
nextflow run main.nf \
  --samplesheet samplesheet.csv \
  --output_path ./results \
  --email you@example.com
```

## Parameters

| parameter | default | description |
|-----------|---------|-------------|
| `--samplesheet` | *(required)* | Path to input samplesheet CSV |
| `--output_path` | `./` | Directory for results and pipeline reports |
| `--downsample` | `true` | Create a downsampled Seurat object alongside the full one |
| `--cluster` | `true` | Run sketch-based clustering |
| `--cluster_full` | `false` | Run full clustering (slower, more memory) |
| `--score_markers` | `false` | Score cell-type marker genes |
| `--email` | *(none)* | Email address for completion/failure notification |

## Profiles

| profile | description |
|---------|-------------|
| `standard` | Local execution, no container |
| `docker` | Local execution with Docker |
| `singularity` | Local execution with Singularity |
| `apptainer` | Local execution with Apptainer |
| `oscer` | SLURM execution on OSCER with Apptainer |

```bash
nextflow run main.nf -profile oscer --samplesheet samplesheet.csv
```

## Outputs

Results are written to `<output_path>/results/<sample>/` and pipeline metadata (trace, timeline, DAG, report) to `<output_path>/pipeline_info/`.