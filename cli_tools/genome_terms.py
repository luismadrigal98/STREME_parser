"""
Shared helpers for genome terminology and reference handling.

The toolkit's canonical term is "genome" (a Penstemon/Mimulus assembly, an
ecotype, an inbred line, ...). Earlier versions used "line" / "Line", and the
analyzers still use those names internally. These helpers let every tool accept
both the canonical column names and the legacy ones, so old Mimulus data keeps
working without migration.

Canonical column names:
  - consolidated sites / motif tables:  genome   (legacy alias: line)
  - long-format expression tables:       Genome   (legacy alias: Line)
"""

DEFAULT_REFERENCE = "IM767"

# Consolidated-sites / motif tables.
_MOTIF_ALIASES = ("genome", "Genome", "line", "Line")
# Long-format expression tables.
_EXPR_ALIASES = ("Genome", "genome", "Line", "line")


def normalize_motif_genome_column(df, target="line"):
    """
    Ensure a consolidated-sites DataFrame exposes `target`.

    Accepts files written with the canonical 'genome' column or the legacy
    'line' column. Default target is the legacy 'line' so the analyzers'
    internal logic stays unchanged.
    """
    if target in df.columns:
        return df
    for alias in _MOTIF_ALIASES:
        if alias != target and alias in df.columns:
            return df.rename(columns={alias: target})
    return df


def normalize_expression_genome_column(df, target="Line"):
    """Ensure a long-format expression DataFrame exposes `target` (default 'Line')."""
    if target in df.columns:
        return df
    for alias in _EXPR_ALIASES:
        if alias != target and alias in df.columns:
            return df.rename(columns={alias: target})
    return df


def row_genome(row):
    """Read the genome value from a csv.DictReader row, canonical or legacy."""
    for key in _MOTIF_ALIASES:
        if key in row and row[key] not in (None, ""):
            return row[key]
    return ""


def add_reference_argument(parser, default=DEFAULT_REFERENCE):
    """
    Add the canonical --reference-genome flag (with --reference-line kept as a
    deprecated alias). Stored on `args.reference_line` so existing analyzer code
    that reads that attribute is unaffected.
    """
    parser.add_argument(
        "--reference-genome", "--reference-line", "-r",
        dest="reference_line", default=default,
        help="Reference genome used as the baseline for relative analysis "
             "(default: %(default)s). '--reference-line' is a deprecated alias.",
    )


def is_wide_expression(df):
    """
    True when the expression table is wide (one column per genome) rather than
    long (Gene / Genome / Expression). Generic: does not assume 'IM' names.
    """
    lower = {c.lower() for c in df.columns}
    if "expression" in lower:
        return False
    return len(df.columns) >= 3


def wide_genome_columns(df, exclude=("LRTadd",)):
    """Genome columns of a wide expression table (everything but gene + excludes)."""
    gene_col = df.columns[0]
    return [c for c in df.columns if c != gene_col and c not in exclude]


def infer_zero_reference(df, genome_cols):
    """
    Return the genome column whose values are all ~0 (the relative baseline),
    or None. Lets relative analysis auto-detect the reference when unspecified.
    """
    import numpy as np
    for col in genome_cols:
        try:
            if np.allclose(df[col].astype(float).values, 0.0, atol=1e-6):
                return col
        except (TypeError, ValueError):
            continue
    return None
