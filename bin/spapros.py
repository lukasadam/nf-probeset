import argparse
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
import yaml

SCRIPT_DIR = str(Path(__file__).resolve().parent)
if SCRIPT_DIR in sys.path:
	sys.path.remove(SCRIPT_DIR)

if not hasattr(np, "in1d"):
	np.in1d = np.isin

import spapros as sp


def parse_args() -> argparse.Namespace:
	parser = argparse.ArgumentParser(
		description="Run a configurable SPAProS probeset selection workflow."
	)
	parser.add_argument(
		"--input",
		type=Path,
		help="Path to the input AnnData file (.h5ad) with raw counts in X.",
	)
	parser.add_argument(
		"--dataset",
		default="pbmc3k_annotated",
		help="Scanpy dataset loader name when --input is not supplied.",
	)
	parser.add_argument(
		"--celltype-target-key",
		default="celltype",
		help="Observation key in adata.obs used as the cell type target.",
	)
	parser.add_argument(
		"--embedding-key",
		default="X_umap",
		help="Embedding key expected in adata.obsm.",
	)
	parser.add_argument(
		"--n-top-genes",
		type=int,
		default=1000,
		help="Number of highly variable genes to compute.",
	)
	parser.add_argument(
		"--hvg-flavor",
		default="cell_ranger",
		help="Flavor passed to scanpy.pp.highly_variable_genes.",
	)
	parser.add_argument(
		"--probeset-size",
		type=int,
		default=20,
		help="Target number of genes for the selected probeset.",
	)
	parser.add_argument(
		"--selector-verbosity",
		type=int,
		default=1,
		help="Verbosity passed to spapros.se.ProbesetSelector.",
	)
	parser.add_argument(
		"--scanpy-verbosity",
		type=int,
		default=0,
		help="scanpy.settings.verbosity value.",
	)
	parser.add_argument(
		"--skip-normalize-total",
		action="store_true",
		help="Skip scanpy.pp.normalize_total.",
	)
	parser.add_argument(
		"--skip-log1p",
		action="store_true",
		help="Skip scanpy.pp.log1p.",
	)
	parser.add_argument(
		"--skip-hvg",
		action="store_true",
		help="Skip scanpy.pp.highly_variable_genes.",
	)
	parser.add_argument(
		"--save-dir",
		type=Path,
		help="Optional directory passed to SPAProS for intermediate selector outputs.",
	)
	parser.add_argument(
		"--selector-config",
		type=Path,
		help=(
			"Path to a YAML or JSON file containing kwargs for "
			"spapros.se.ProbesetSelector."
		),
	)
	parser.add_argument(
		"--selector-kwargs",
		default="{}",
		help=(
			"JSON object merged into spapros.se.ProbesetSelector kwargs. "
			"Values here override the wrapper defaults."
		),
	)
	parser.add_argument(
		"--output-dir",
		type=Path,
		default=Path("."),
		help="Directory for generated outputs.",
	)
	parser.add_argument(
		"--probeset-output",
		default="probeset.csv",
		help="Output CSV filename for the selected probeset.",
	)
	parser.add_argument(
		"--dotplot-output",
		default="dotplot.pdf",
		help="Output filename for the masked dotplot.",
	)
	parser.add_argument(
		"--evaluation-summary-output",
		default="evaluation_summary.csv",
		help="Output filename for the evaluation summary CSV.",
	)
	parser.add_argument(
		"--evaluation-summary-plot",
		default="evaluation_summary.pdf",
		help="Output filename for the evaluation summary plot.",
	)
	return parser.parse_args()


def parse_selector_kwargs(raw_value: str) -> dict:
	try:
		selector_kwargs = json.loads(raw_value)
	except json.JSONDecodeError as exc:
		raise ValueError("--selector-kwargs must be valid JSON.") from exc
	if not isinstance(selector_kwargs, dict):
		raise ValueError("--selector-kwargs must decode to a JSON object.")
	return selector_kwargs


def load_selector_config(config_path: Path | None) -> dict:
	if config_path is None:
		return {}
	with config_path.open("r", encoding="utf-8") as handle:
		if config_path.suffix.lower() == ".json":
			selector_kwargs = json.load(handle)
		else:
			selector_kwargs = yaml.safe_load(handle)
	if selector_kwargs is None:
		return {}
	if not isinstance(selector_kwargs, dict):
		raise ValueError("--selector-config must contain a mapping/object.")
	return selector_kwargs


def load_pbmc3k_annotated():
	adata = sc.datasets.pbmc3k()
	adata_processed = sc.datasets.pbmc3k_processed()
	adata = adata[adata_processed.obs_names, adata_processed.var_names].copy()
	adata.obs["celltype"] = adata_processed.obs["louvain"]
	adata.obsm["X_umap"] = adata_processed.obsm["X_umap"]
	return adata


def load_adata(input_path: Path | None, dataset_name: str):
	if input_path is not None:
		return sc.read_h5ad(input_path)
	if dataset_name == "pbmc3k_annotated":
		return load_pbmc3k_annotated()
	try:
		dataset_loader = getattr(sc.datasets, dataset_name)
	except AttributeError as exc:
		raise ValueError(f"Unknown scanpy dataset: {dataset_name}") from exc
	return dataset_loader()


def validate_adata(adata, celltype_key: str, embedding_key: str) -> None:
	if adata.n_obs == 0 or adata.n_vars == 0:
		raise ValueError("Input AnnData is empty.")
	if adata.X is None:
		raise ValueError("Input AnnData must contain raw counts in X.")
	if celltype_key not in adata.obs:
		raise KeyError(
			f"Cell type key '{celltype_key}' was not found in adata.obs."
		)
	if embedding_key not in adata.obsm:
		raise KeyError(
			f"Embedding key '{embedding_key}' was not found in adata.obsm."
		)


def main() -> None:
	args = parse_args()
	sc.settings.verbosity = args.scanpy_verbosity

	output_dir = args.output_dir
	output_dir.mkdir(parents=True, exist_ok=True)
	sc.logging.print_header()
	print(f"spapros=={getattr(sp, '__version__', 'unknown')}")

	adata = load_adata(args.input, args.dataset)
	validate_adata(adata, args.celltype_target_key, args.embedding_key)

	if not args.skip_normalize_total:
		sc.pp.normalize_total(adata)
	if not args.skip_log1p:
		sc.pp.log1p(adata)
	if not args.skip_hvg:
		sc.pp.highly_variable_genes(
			adata,
			flavor=args.hvg_flavor,
			n_top_genes=args.n_top_genes,
		)

	save_dir = str(args.save_dir) if args.save_dir is not None else None
	selector_kwargs = {
		"n": args.probeset_size,
		"celltype_key": args.celltype_target_key,
		"verbosity": args.selector_verbosity,
		"save_dir": save_dir,
	}
	selector_kwargs.update(load_selector_config(args.selector_config))
	selector_kwargs.update(parse_selector_kwargs(args.selector_kwargs))

	# Run probeset selection
	selector = sp.se.ProbesetSelector(adata, **selector_kwargs)
	selector.select_probeset()

	# Save selected probeset to CSV
	selector.probeset.to_csv(output_dir / args.probeset_output, index=False)

	# Generate and save masked dotplot
	sp.pl.masked_dotplot(
		adata,
		selector,
		ct_key=args.celltype_target_key,
		save=str(output_dir / args.dotplot_output),
	)

	# Infer reference probesets
	reference_sets = sp.se.select_reference_probesets(adata, n=args.probeset_size)

	# Evaluate the selected probeset against the inferred references
	evaluator = sp.ev.ProbesetEvaluator(adata, verbosity=2, results_dir=None)

	# Get the selected probeset genes and evaluate against references (DE, PCA, random)
	probeset = list(selector.probeset.index[selector.probeset.selection])
	evaluator.evaluate_probeset(probeset, set_id="nf-probeset")

	# Evaluate reference sets
	for set_id, df in reference_sets.items():
		gene_set = df[df["selection"]].index.to_list()
		evaluator.evaluate_probeset(gene_set, set_id=set_id)
	
	# Save evaluation summary
	evaluator.summary_results.to_csv(output_dir / args.evaluation_summary_output, index=False)

	# Generate evaluation summary plot
	evaluator.plot_summary(save=output_dir / args.evaluation_summary_plot)

if __name__ == "__main__":
	main()