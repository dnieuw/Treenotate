# Treenotate

Treenotate annotates Newick trees using tip metadata and can optionally add APOBEC3-style mutation branch labels from IQ-TREE ancestral-state output.

## Installation (Conda/Mamba)

Use `mamba` (or `conda`) to create an isolated environment, then install the tool.

### Recommended: create from environment file

From the repository directory:

```bash
mamba env create -f environment.yml
mamba activate treenotate
```

### 2. Install Treenotate

From this repository directory:

```bash
pip install .
```

This installs the command-line tool `treenotate`.

### 3. Verify installation

```bash
treenotate --version
treenotate --help
```

## Usage

Treenotate supports two modes:

1. Metadata-only annotation (tip annotations only)
2. Metadata + APOBEC branch annotation (requires both alignment and state file)

### Required inputs

- `--tree-file`: input Newick tree
- `--metadata-file`: tabular metadata file for tip annotations
- `--output-tree`: output annotated Newick tree path

### Metadata options

- `--metadata-delimiter` (default `,`)
- `--tip-column` (default `old_tiplabel`)

### Optional APOBEC analysis

To enable APOBEC branch labels, provide both:

- `--alignment`: FASTA alignment used to make the tree
- `--state-file`: IQ-TREE `.state` file (with `Site`, `Node`, `State` columns)

If only one of these is provided, Treenotate exits with an error.

## Examples

### Metadata-only mode

```bash
treenotate \
	--tree-file data/tree.treefile \
	--metadata-file data/metadata.csv \
	--tip-column old_tiplabel \
	--output-tree results/tree.annotated.treefile
```

### Metadata + APOBEC mode

```bash
treenotate \
	--tree-file data/tree.treefile \
	--metadata-file data/metadata.csv \
	--alignment data/alignment.fasta \
	--state-file data/tree.state \
	--output-tree results/tree.apobec.annotated.treefile
```

## Input format notes

- Tip names in the tree must match values in the metadata tip column.
- Metadata file is read as text; values are cleaned and written into BEAST-style comment blocks.
- IQTREE state file can be generated with IQTREE by using the `-asr` option to 'Write ancestral sequences (by empirical Bayesian method) for all nodes of the tree to .state file.'
