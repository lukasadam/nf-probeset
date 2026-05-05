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

Results are written to `results/` by default.

## Useful Parameters

- `--dataset`
- `--input`
- `--probeset-size`
- `--output-dir`
- `--container`