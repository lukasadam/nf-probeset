# nf-probeset

Minimal Nextflow pipeline for generating a probe set with `spapros`.

## Requirements

- Nextflow
- Python environment at `.venv` with packages from `requirements.txt`, or a container image passed via `--container`

## Run

Run with the default bundled dataset:

```bash
nextflow run main.nf
```

Run with an input AnnData file:

```bash
nextflow run main.nf --input path/to/data.h5ad
```

Prefer putting `spapros.se.ProbesetSelector(...)` kwargs in a YAML config file:

```bash
nextflow run main.nf --selector_config selector.example.yml
```

The config file is loaded and passed through to `ProbesetSelector`, overriding the wrapper defaults where needed.

For quick one-off overrides, you can still pass an inline map:

```bash
nextflow run main.nf --selector_kwargs '{"n_pca_genes":250,"forest_hparams":{"n_trees":100}}'
```

Results are written to `results/` by default.

## Container Usage

The GitHub Actions workflow publishes a Docker image to GHCR. It does not build a `.sif` artifact in CI.

If you need a Singularity/Apptainer image, pull it from GHCR on a machine that has Apptainer installed:

```bash
apptainer pull docker://ghcr.io/lukasadam/nf-probeset:latest
```

You can also point Nextflow at the published container image:

```bash
nextflow run main.nf -profile singularity --container docker://ghcr.io/lukasadam/nf-probeset:latest
```

## Useful Parameters

- `--dataset`
- `--input`
- `--probeset-size`
- `--output-dir`
- `--container`