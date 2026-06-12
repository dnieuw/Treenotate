#!/usr/bin/env python3

import argparse
import logging
import re
import os
from collections import Counter
from pathlib import Path

import pandas as pd
import numpy as np
from Bio import Phylo, SeqIO


LOGGER = logging.getLogger("Treenotate")
_NUM_RE = re.compile(r"^-?\d+(\.\d+)?([eE]-?\d+)?$")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Treenotate: annotate a phylogenetic tree with tip metadata and optional APOBEC3 branch labels."
    )
    parser.add_argument(
        "-a",
        "--alignment",
        required=False,
        help="Input FASTA alignment used to provide tip states.",
    )
    parser.add_argument(
        "-t",
        "--tree-file",
        required=True,
        help="Input Newick tree file to annotate.",
    )
    parser.add_argument(
        "-s",
        "--state-file",
        required=False,
        help="Input IQTREE state file containing Site/Node/State columns. (use -asr option in IQTREE to 'Write ancestral sequences (by empirical Bayesian method) for all nodes of the tree to .state file.')",
    )
    parser.add_argument(
        "-m",
        "--metadata-file",
        required=True,
        help="Input metadata table for tree tip annotations.",
    )
    parser.add_argument(
        "-o",
        "--output-tree",
        required=True,
        help="Output path for annotated Newick tree.",
    )
    parser.add_argument(
        "-d",
        "--metadata-delimiter",
        default=",",
        help="Delimiter for metadata file (default: ',').",
    )
    parser.add_argument(
        "-c",
        "--tip-column",
        default="old_tiplabel",
        help="Metadata column containing tip names matching tree terminals.",
    )
    parser.add_argument(
        "-l",
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        help="Logging verbosity (default: INFO).",
    )
    parser.add_argument(
        "-L",
        "--log-file",
        default=None,
        help="Optional path to write logs in addition to stderr.",
    )
    parser.add_argument(
        "--version",
        action="version",
        version="Treenotate 1.0.0",
        help="Show program version and exit.",
    )
    return parser.parse_args()


def configure_logging(log_level: str, log_file: str | None) -> None:
    level = getattr(logging, log_level.upper(), logging.INFO)
    handlers: list[logging.Handler] = [logging.StreamHandler()]
    if log_file:
        handlers.append(logging.FileHandler(log_file))

    logging.basicConfig(
        level=level,
        format="%(asctime)s | %(levelname)s | %(name)s | %(message)s",
        handlers=handlers,
    )


def load_annotations_df(
    annotation_path: str | Path,
    *,
    delimiter: str,
    tip_column: str,
) -> pd.DataFrame:
    """Read and clean annotation file into a DataFrame indexed by tip id."""
    annotation_path = Path(annotation_path)
    if not annotation_path.exists() or not annotation_path.is_file():
        raise FileNotFoundError(f"Could not find annotation file: {annotation_path}")

    df = pd.read_csv(annotation_path, sep=delimiter, dtype=str).fillna("")
    if tip_column not in df.columns:
        raise ValueError(f"tip_column={tip_column!r} not in headers: {list(df.columns)}")

    for col in df.columns:
        df[col] = df[col].astype(str).str.strip()

    if df[tip_column].eq("").any():
        raise ValueError(f"Missing tips in column {tip_column!r}")

    return df.set_index(tip_column, verify_integrity=True)


def beast_fmt_value(raw: str) -> str:
    """Format a value for `[&key=value]` blocks."""
    s = raw.strip()
    if _NUM_RE.match(s):
        return s
    if s in ("true", "false"):
        return s
    if s.startswith("{") and s.endswith("}"):
        return s
    s = s.replace("\\", "\\\\").replace('"', '\\"')
    return f'"{s}"'


def build_beast_block_from_row(row: pd.Series) -> str:
    """Build `key=value,...` from one annotation DataFrame row."""
    parts = []
    for key, value in row.items():
        if value is None:
            continue
        vv = str(value).strip()
        if vv == "":
            continue
        parts.append(f"{key}={beast_fmt_value(vv)}")
    return ",".join(parts)
    

def is_variable(row: np.array) -> bool:
    """A site is variable if more than one A/C/G/T state is present."""
    unique = np.unique(row)
    unique = [k for k in unique if k in ["A","C","G","T"]]
    return len(unique) > 1


def create_dimers(states_array: np.ndarray, site: int) -> np.ndarray:
    """Create dimers at positions site and site+1 for all nodes."""
    return np.char.add(states_array[site,:], states_array[site+1,:])


def check_apobec_mutations(
    node1: str,
    node2: str,
    states_df_variable: pd.DataFrame,
    dimers_df: pd.DataFrame,
    dimers_prev_df: pd.DataFrame,
) -> Counter:
    """Count APOBEC-consistent vs non-APOBEC branch changes between two nodes."""
    branch = states_df_variable[[node1, node2]]
    variable_positions = branch.index[branch.apply(is_variable, axis=1, raw=True)]

    branch_dimers = pd.concat(
        [
            dimers_df[[node1, node2]].loc[variable_positions].rename(columns={node1: "x", node2: "y"}),
            dimers_prev_df[[node1, node2]].loc[variable_positions].rename(
                columns={node1: "i", node2: "j"}
            ),
        ],
        axis=1,
    )

    variable_dimers = branch_dimers[
        (branch_dimers["x"] != branch_dimers["y"]) | (branch_dimers["i"] != branch_dimers["j"])
    ]

    mask_apobec = (
        ((variable_dimers["x"] == "TC") & (variable_dimers["y"] == "TT"))
        | ((variable_dimers["x"] == "GA") & (variable_dimers["y"] == "AA"))
        | ((variable_dimers["i"] == "TC") & (variable_dimers["j"] == "TT"))
        | ((variable_dimers["i"] == "GA") & (variable_dimers["j"] == "AA"))
    )
    branch_apobec3 = variable_dimers[mask_apobec]

    return Counter({"APOBEC3": len(branch_apobec3), "non-APOBEC3": len(variable_positions) - len(branch_apobec3)})

def annotate_clade(
    clade: Phylo.Newick.Clade,
    ann_df: pd.DataFrame,
    states_df_variable: pd.DataFrame | None = None,
    dimers_df: pd.DataFrame | None = None,
    dimers_prev_df: pd.DataFrame | None = None,
) -> None:
    """Recursively annotate tree clades with metadata and APOBEC branch labels."""
    strip_digits = lambda s: re.sub(r"/\d+$", "", s)
    strip_node = lambda s: re.sub(r"^Node\d+/", "", s)

    if clade.is_terminal():
        annotations_row = ann_df.loc[clade.name]
        annotations = build_beast_block_from_row(annotations_row)
        if clade.comment == "" or clade.comment is None:
            clade.comment = f"&{annotations}"
        else:
            clade.comment = f"{clade.comment},{annotations}"
        return

    apobec_enabled = (
        states_df_variable is not None and dimers_df is not None and dimers_prev_df is not None
    )

    current_node = strip_digits(clade.name)
    clade.name = strip_node(clade.name)

    for child in clade.clades:
        if apobec_enabled:
            child_node = strip_digits(child.name)
            counts = check_apobec_mutations(
                current_node, child_node, states_df_variable, dimers_df, dimers_prev_df
            )
            if counts.total() > 0:
                child.comment = f'&APOBEC3_label="{counts["APOBEC3"]}/{counts.total()}"'
        annotate_clade(child, ann_df, states_df_variable, dimers_df, dimers_prev_df)

def read_states_array(state_file: str) -> np.ndarray:
    LOGGER.info("Reading state file: %s", state_file)

    #Find last state to determine array shape
    with open(state_file, 'rb') as f:
        try:
            f.seek(-2, os.SEEK_END)
            while f.read(1) != b'\n':
                f.seek(-2, os.SEEK_CUR)
        except OSError:
            f.seek(0)
        last_line = f.readline().decode()
        n_states = int(last_line.split('\t')[1])

    array_list = []
    node_names = []
    with open(state_file, 'r') as f:
        current_node = None
        for line in f:
            if line.startswith("#"):
                continue
            if line.startswith("Node\t"):
                continue

            try:
                node,site,state,probA,probC,probG,probT = line.rstrip("\n").split("\t")
            except:
                raise ValueError(f"Malformed line in state file: {line.strip()}")

            if not current_node:
                current_node = node
                node_names.append(node)
                states_array = np.empty(n_states, dtype='U1')
            elif node != current_node:
                array_list.append(states_array)

                current_node = node
                node_names.append(node)
                states_array = np.empty(n_states, dtype='U1')
            else:
                states_array[int(site)-1] = state
        #Add last array
        array_list.append(states_array)

    states_array = np.stack(array_list, axis=-1)

    LOGGER.info("Loaded %d sites from %d node columns", states_array.shape[0], states_array.shape[1])
    return states_array, node_names

def read_fasta_as_array(alignment_path: str) -> np.ndarray:
    LOGGER.info("Processing FASTA alignment: %s", alignment_path)
    fasta_columns = []
    fasta_names = []
    for record in SeqIO.parse(alignment_path, "fasta"):
        seq_column = np.array(list(record.seq))
        fasta_columns.append(seq_column)
        fasta_names.append(record.id)

    if not fasta_columns:
        raise ValueError(f"No sequences found in FASTA alignment: {alignment_path}")

    fasta_array = np.stack(fasta_columns, axis=1)
    LOGGER.info("Loaded %d sites from %d tip sequences", fasta_array.shape[0], fasta_array.shape[1])
    return fasta_array, fasta_names

def generate_states_array(alignment_path: str, state_file: str) -> np.ndarray:
    states_array, node_names = read_states_array(state_file)
    fasta_array, fasta_names = read_fasta_as_array(alignment_path)
    
    if fasta_array.shape[0] != states_array.shape[0]:
        raise ValueError(
            f"Number of sites in FASTA ({fasta_array.shape[0]}) does not match state file ({states_array.shape[0]})."
        )

    combined_array = np.concatenate([fasta_array, states_array], axis=1)
    combined_names = fasta_names + node_names
    LOGGER.info("Combined state matrix shape: %s", combined_array.shape)
    return combined_array, combined_names

def run_workflow(args: argparse.Namespace) -> None:
    has_alignment = bool(args.alignment)
    has_state_file = bool(args.state_file)
    if has_alignment != has_state_file:
        raise ValueError(
            "To enable APOBEC annotation, provide both --alignment and --state-file. "
            "For metadata-only annotation, provide neither."
        )

    apobec_enabled = has_alignment and has_state_file

    states_df_variable: pd.DataFrame | None = None
    dimers_df: pd.DataFrame | None = None
    dimers_prev_df: pd.DataFrame | None = None
    if apobec_enabled:
        LOGGER.info("APOBEC mode enabled: reading state and alignment inputs")
        states_array, states_names = generate_states_array(args.alignment, args.state_file)

        LOGGER.info("Finding variable sites")
        variable_positions = np.where(np.apply_along_axis(is_variable, 1, states_array))[0]
        states_df_variable = pd.DataFrame(states_array[variable_positions,:], columns=states_names)
        LOGGER.info("Identified %d variable positions", len(variable_positions))

        LOGGER.info("Creating dimers")
        dimers_df = pd.DataFrame(np.stack([create_dimers(states_array, i) for i in variable_positions], axis=1).T, columns=states_names)
        dimers_prev_df = pd.DataFrame(np.stack([create_dimers(states_array, i) for i in variable_positions-1], axis=1).T, columns=states_names)
    else:
        LOGGER.info("Metadata-only mode enabled: APOBEC branch annotation disabled")

    LOGGER.info("Reading tree: %s", args.tree_file)
    tree = Phylo.read(args.tree_file, "newick")

    LOGGER.info("Reading metadata: %s", args.metadata_file)
    ann_df = load_annotations_df(
        args.metadata_file,
        delimiter=args.metadata_delimiter,
        tip_column=args.tip_column,
    )

    LOGGER.info("Annotating tree")
    annotate_clade(tree.clade, ann_df, states_df_variable, dimers_df, dimers_prev_df)

    tree.clade.name = ""
    output_path = Path(args.output_tree)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    LOGGER.info("Writing annotated tree: %s", output_path)
    Phylo.write(tree, str(output_path), "newick", format_branch_length="%1.10f")


def main() -> int:
    args = parse_args()
    configure_logging(args.log_level, args.log_file)
    LOGGER.info("Starting Treenotate workflow")

    try:
        run_workflow(args)
    except Exception:
        LOGGER.exception("Treenotate workflow failed")
        return 1

    LOGGER.info("Treenotate workflow completed successfully")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())