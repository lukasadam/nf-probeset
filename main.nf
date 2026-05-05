nextflow.enable.dsl=2

params.input = null
params.dataset = 'pbmc3k_annotated'
params.celltype_target_key = 'celltype'
params.embedding_key = 'X_umap'
params.n_top_genes = 1000
params.hvg_flavor = 'cell_ranger'
params.probeset_size = 20
params.selector_verbosity = 1
params.scanpy_verbosity = 0
params.skip_normalize_total = false
params.skip_log1p = false
params.skip_hvg = false
params.save_dir = null
params.selector_config = null
params.selector_kwargs = [:]
params.output_dir = 'results'
params.probeset_output = 'probeset.csv'
params.dotplot_output = 'dotplot.pdf'
params.evaluation_summary_output = 'evaluation_summary.csv'
params.evaluation_summary_plot = 'evaluation_summary.pdf'
params.python_command = "${projectDir}/.venv/bin/python"
params.container = null

process RUN_SPAPROS {
    tag "spapros"
    publishDir params.output_dir, mode: 'copy'

    output:
    path params.probeset_output
    path params.dotplot_output
    path params.evaluation_summary_plot
    path params.evaluation_summary_output

    script:
    def cli = []
    def pythonCommand = params.container ? 'python' : params.python_command
    def scriptPath = params.container ? '/opt/nf-probeset/bin/spapros.py' : "${projectDir}/bin/spapros.py"
    cli << pythonCommand
    cli << scriptPath
    if( params.input ) cli += ['--input', params.input]
    else cli += ['--dataset', params.dataset]
    cli += ['--celltype-target-key', params.celltype_target_key]
    cli += ['--embedding-key', params.embedding_key]
    cli += ['--n-top-genes', params.n_top_genes.toString()]
    cli += ['--hvg-flavor', params.hvg_flavor]
    cli += ['--probeset-size', params.probeset_size.toString()]
    cli += ['--selector-verbosity', params.selector_verbosity.toString()]
    cli += ['--scanpy-verbosity', params.scanpy_verbosity.toString()]
    if( params.selector_config ) cli += ['--selector-config', params.selector_config]
    cli += ['--selector-kwargs', groovy.json.JsonOutput.toJson(params.selector_kwargs)]
    cli += ['--output-dir', '.']
    cli += ['--probeset-output', params.probeset_output]
    cli += ['--dotplot-output', params.dotplot_output]
    if( params.save_dir ) cli += ['--save-dir', params.save_dir]
    if( params.skip_normalize_total ) cli << '--skip-normalize-total'
    if( params.skip_log1p ) cli << '--skip-log1p'
    if( params.skip_hvg ) cli << '--skip-hvg'

    """
    ${cli.collect { "'${it.toString().replace("'", "'\\''")}'" }.join(' ')}
    """
}

workflow {
    RUN_SPAPROS()
}
