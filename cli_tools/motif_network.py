#!/usr/bin/env python3
"""
Motif co-occurrence networks and gene clustering by regulatory profile.

Reads a consolidated STREME-sites TSV (from streme_sites_consolidator) and
produces two complementary explorations of the data — no expression data
required:

  1. Motif co-occurrence network: pairs of motifs scored by Jaccard / lift /
     Fisher one-sided p (with BH-corrected q-values); thresholded into a
     network whose connected components are interpreted as "modules" of
     motifs that tend to appear together in the same promoters. Exported as
     GraphML for direct import into Cytoscape / Gephi.
  2. Gene clustering by motif profile: hierarchical clustering of the gene ×
     motif binary presence matrix (Jaccard distance, average linkage), with
     per-cluster motif fingerprints (Fisher one-sided enrichment vs the
     genome background, BH-corrected within cluster).

A light positional summary (start_pos / relative_position_fraction
distribution per motif) is also produced. After TOMTOM annotation, the
graphml node IDs are STREME motif IDs which can be joined to TF names via
the tomtom.tsv output for biological interpretation.

Outputs (in --output directory):
  motif_cooccurrence_edges.tsv         qualifying motif pairs with stats
  motif_cooccurrence_network.graphml   importable into Cytoscape / Gephi
  motif_modules.tsv                    motif_id -> module_id (connected comps)
  gene_clusters.tsv                    gene_id -> cluster_id (NA for unclustered)
  gene_cluster_fingerprints.tsv        per cluster, enriched motifs
  motif_positional_summary.tsv         per motif, position-from-promoter-start stats
  network_summary.txt                  human-readable run summary
"""

import argparse
import os
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import pdist

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import genome_terms  # noqa: E402  (after sys.path adjustment)


# --------------------------------------------------------------------------- #
# I/O
# --------------------------------------------------------------------------- #

def load_consolidated(motif_file):
    df = pd.read_csv(motif_file, sep="\t")
    df = genome_terms.normalize_motif_genome_column(df)
    required = {"consolidated_motif_id", "gene_id"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(
            f"Consolidated file missing required columns: {sorted(missing)}"
        )
    return df


def build_presence_matrix(df, min_motif_sites=10):
    """gene × motif binary DataFrame, restricted to motifs present in >= N genes."""
    motif_gene_counts = df.groupby("consolidated_motif_id")["gene_id"].nunique()
    keep = motif_gene_counts[motif_gene_counts >= min_motif_sites].index
    sub = df[df["consolidated_motif_id"].isin(keep)]
    presence = (
        sub.groupby(["gene_id", "consolidated_motif_id"]).size().unstack(fill_value=0)
    )
    return (presence > 0).astype(np.int8)


# --------------------------------------------------------------------------- #
# Multiple-testing helper (BH q-values, monotone)
# --------------------------------------------------------------------------- #

def bh_qvalues(pvals):
    """Benjamini–Hochberg adjusted q-values, monotone (Series indexed like input)."""
    p = pd.Series(pvals).astype(float)
    n = len(p)
    if n == 0:
        return p
    order = p.sort_values(kind="mergesort").index
    ranks = np.arange(1, n + 1)
    sorted_p = p.loc[order].values
    raw = np.minimum(sorted_p * n / ranks, 1.0)
    # Monotone non-decreasing from the right: cummin of the reversed array.
    monotone = np.minimum.accumulate(raw[::-1])[::-1]
    out = pd.Series(monotone, index=order)
    return out.reindex(p.index)


# --------------------------------------------------------------------------- #
# Motif co-occurrence
# --------------------------------------------------------------------------- #

def motif_cooccurrence(presence, min_jaccard=0.1, min_lift=2.0, fdr=0.05):
    """
    All-pairs motif co-occurrence with Jaccard, lift, Fisher one-sided p
    (enrichment), and BH q-values.

    Returns:
        all_edges_df  — every pair with non-zero co-occurrence
        kept_edges_df — pairs passing (jaccard>=min_jaccard OR lift>=min_lift) AND q<=fdr,
                        sorted by jaccard desc.
    """
    motifs = list(presence.columns)
    M = presence.values  # int8 (n_genes, n_motifs)
    n_genes = M.shape[0]
    per_motif_counts = M.sum(axis=0).astype(np.int64)

    # Vectorised pairwise co-occurrence: A^T A gives the n_both for every pair.
    both = (M.astype(np.int32).T @ M.astype(np.int32))

    rows = []
    for i in range(len(motifs)):
        ci = int(per_motif_counts[i])
        if ci == 0:
            continue
        for j in range(i + 1, len(motifs)):
            cj = int(per_motif_counts[j])
            if cj == 0:
                continue
            b = int(both[i, j])
            if b == 0:
                continue
            union = ci + cj - b
            jaccard = b / union if union > 0 else 0.0
            lift = (b * n_genes) / (ci * cj)
            # Fisher 2x2: [[both, ci-both], [cj-both, n-ci-cj+both]], one-sided greater
            a = b
            bb = ci - b
            c = cj - b
            d = n_genes - ci - cj + b
            try:
                _, p = stats.fisher_exact([[a, bb], [c, d]], alternative="greater")
            except Exception:
                p = 1.0
            rows.append({
                "motif_a": motifs[i],
                "motif_b": motifs[j],
                "count_a": ci,
                "count_b": cj,
                "count_both": b,
                "n_genes": n_genes,
                "jaccard": jaccard,
                "lift": lift,
                "fisher_p": p,
            })

    if not rows:
        cols = ["motif_a", "motif_b", "count_a", "count_b", "count_both",
                "n_genes", "jaccard", "lift", "fisher_p", "fisher_q"]
        return pd.DataFrame(columns=cols), pd.DataFrame(columns=cols)

    all_edges = pd.DataFrame(rows)
    all_edges["fisher_q"] = bh_qvalues(all_edges["fisher_p"])
    keep = ((all_edges["jaccard"] >= min_jaccard) | (all_edges["lift"] >= min_lift)) \
        & (all_edges["fisher_q"] <= fdr)
    kept = all_edges[keep].sort_values(["jaccard", "lift"], ascending=False).reset_index(drop=True)
    return all_edges, kept


def connected_modules(edges_df, all_motifs):
    """Union-find connected components on the thresholded edge set.

    Returns dict motif -> module_id. Singletons get module_id 0; non-trivial
    modules are numbered 1..K in descending size order.
    """
    parent = {m: m for m in all_motifs}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(x, y):
        rx, ry = find(x), find(y)
        if rx != ry:
            parent[ry] = rx

    for _, e in edges_df.iterrows():
        if e["motif_a"] in parent and e["motif_b"] in parent:
            union(e["motif_a"], e["motif_b"])

    groups = defaultdict(list)
    for m in all_motifs:
        groups[find(m)].append(m)

    ordered = sorted(groups.values(), key=len, reverse=True)
    out = {}
    next_id = 1
    for members in ordered:
        if len(members) == 1:
            out[members[0]] = 0
        else:
            for m in members:
                out[m] = next_id
            next_id += 1
    return out


def write_graphml(edges_df, nodes_meta, path):
    """Minimal GraphML writer (no networkx dependency)."""
    def esc(s):
        return str(s).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;").replace('"', "&quot;")

    lines = [
        '<?xml version="1.0" encoding="UTF-8"?>',
        '<graphml xmlns="http://graphml.graphdrawing.org/xmlns">',
        '  <key id="count" for="node" attr.name="count" attr.type="int"/>',
        '  <key id="module" for="node" attr.name="module" attr.type="int"/>',
        '  <key id="jaccard" for="edge" attr.name="jaccard" attr.type="double"/>',
        '  <key id="lift" for="edge" attr.name="lift" attr.type="double"/>',
        '  <key id="qvalue" for="edge" attr.name="qvalue" attr.type="double"/>',
        '  <graph edgedefault="undirected">',
    ]
    for node, (count, module) in nodes_meta.items():
        lines.append(
            f'    <node id="{esc(node)}">'
            f'<data key="count">{int(count)}</data>'
            f'<data key="module">{int(module)}</data></node>'
        )
    for _, e in edges_df.iterrows():
        lines.append(
            f'    <edge source="{esc(e["motif_a"])}" target="{esc(e["motif_b"])}">'
            f'<data key="jaccard">{e["jaccard"]:.4f}</data>'
            f'<data key="lift">{e["lift"]:.3f}</data>'
            f'<data key="qvalue">{e["fisher_q"]:.3g}</data></edge>'
        )
    lines.append("  </graph>")
    lines.append("</graphml>")
    with open(path, "w") as fh:
        fh.write("\n".join(lines))


# --------------------------------------------------------------------------- #
# Gene clustering
# --------------------------------------------------------------------------- #

def cluster_genes(presence, min_motifs_per_gene=2, max_genes=5000,
                  n_clusters=None, dist_threshold=0.7, random_state=42):
    """
    Hierarchical clustering of genes on Jaccard distance / average linkage.

    Returns (cluster_df, was_subsampled).
        cluster_df  — every gene in `presence`, with cluster_id (NA when
                      filtered out or unclustered).
        was_subsampled — True if the eligible gene set exceeded max_genes
                         and was randomly subsampled.
    """
    gene_motif_counts = presence.sum(axis=1)
    eligible = gene_motif_counts[gene_motif_counts >= min_motifs_per_gene].index
    sub = presence.loc[eligible]
    was_subsampled = False

    if len(sub) > max_genes:
        rng = np.random.default_rng(random_state)
        chosen = rng.choice(len(sub), size=max_genes, replace=False)
        sub = sub.iloc[chosen]
        was_subsampled = True

    cluster_df = pd.DataFrame({"gene_id": presence.index, "cluster_id": pd.NA})
    cluster_df = cluster_df.set_index("gene_id")

    if len(sub) < 3:
        return cluster_df.reset_index(), was_subsampled

    dists = pdist(sub.values.astype(np.uint8), metric="jaccard")
    Z = linkage(dists, method="average")
    if n_clusters:
        labels = fcluster(Z, t=int(n_clusters), criterion="maxclust")
    else:
        labels = fcluster(Z, t=float(dist_threshold), criterion="distance")

    for gene, lbl in zip(sub.index, labels):
        cluster_df.loc[gene, "cluster_id"] = int(lbl)
    return cluster_df.reset_index(), was_subsampled


def cluster_fingerprints(presence, cluster_df, fdr=0.05, top_n=15, min_cluster_size=3):
    """For each cluster, motifs enriched vs the genome background (Fisher + BH)."""
    valid = cluster_df.dropna(subset=["cluster_id"])
    if valid.empty:
        return pd.DataFrame()

    bg_total = len(presence)
    bg_counts = presence.sum(axis=0)

    out = []
    for cluster_id, group in valid.groupby("cluster_id"):
        cluster_genes = list(set(group["gene_id"]) & set(presence.index))
        n_cluster = len(cluster_genes)
        if n_cluster < min_cluster_size:
            continue
        cluster_counts = presence.loc[cluster_genes].sum(axis=0)

        rows = []
        for motif in presence.columns:
            cc = int(cluster_counts[motif])
            if cc == 0:
                continue
            bc = int(bg_counts[motif])
            a, b = cc, n_cluster - cc
            c = bc - cc
            d = bg_total - n_cluster - c
            if c < 0 or d < 0:
                continue
            try:
                _, p = stats.fisher_exact([[a, b], [c, d]], alternative="greater")
            except Exception:
                p = 1.0
            cluster_freq = cc / n_cluster
            bg_freq = bc / bg_total
            enrichment = cluster_freq / bg_freq if bg_freq > 0 else float("inf")
            rows.append({
                "cluster_id": int(cluster_id),
                "motif_id": motif,
                "cluster_count": cc,
                "cluster_size": n_cluster,
                "cluster_freq": cluster_freq,
                "bg_count": bc,
                "bg_size": bg_total,
                "bg_freq": bg_freq,
                "enrichment": enrichment,
                "fisher_p": p,
            })
        if not rows:
            continue
        cl = pd.DataFrame(rows)
        cl["fisher_q"] = bh_qvalues(cl["fisher_p"])
        top = (cl[cl["fisher_q"] <= fdr]
               .sort_values(["enrichment", "fisher_q"], ascending=[False, True])
               .head(top_n))
        if not top.empty:
            out.append(top)

    return pd.concat(out, ignore_index=True) if out else pd.DataFrame()


# --------------------------------------------------------------------------- #
# Positional summary
# --------------------------------------------------------------------------- #

def positional_summary(df, motifs):
    """Per motif: distribution stats over start_pos and relative_position_fraction."""
    if "start_pos" not in df.columns:
        return pd.DataFrame()
    sub = df[df["consolidated_motif_id"].isin(motifs)]
    if sub.empty:
        return pd.DataFrame()
    pos = (sub.groupby("consolidated_motif_id")["start_pos"]
           .agg(["count", "mean", "median", "std", "min", "max"])
           .reset_index())
    if "relative_position_fraction" in sub.columns:
        rel = (sub.groupby("consolidated_motif_id")["relative_position_fraction"]
               .agg(["mean", "median"])
               .rename(columns={"mean": "rel_pos_mean", "median": "rel_pos_median"})
               .reset_index())
        pos = pos.merge(rel, on="consolidated_motif_id", how="left")
    return pos.sort_values("count", ascending=False)


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #

def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Motif co-occurrence networks + gene clustering on a "
                    "consolidated STREME-sites TSV (no expression data needed).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Defaults: drop motifs with <10 sites, edges need Jaccard>=0.1 or lift>=2,
  # q<=0.05; cluster genes with >=2 motifs at distance 0.7.
  %(prog)s outputs/consolidated_streme_sites.tsv --output network_results/

  # Tighter network, more clusters
  %(prog)s outputs/consolidated_streme_sites.tsv -o net/ \\
      --min-jaccard 0.2 --min-lift 3 --fdr 0.01 --n-clusters 30
""",
    )
    parser.add_argument("consolidated_file",
                        help="Consolidated STREME sites TSV "
                             "(from streme_sites_consolidator)")
    parser.add_argument("--output", "-o", required=True, help="Output directory")
    parser.add_argument("--min-motif-sites", type=int, default=10,
                        help="Drop motifs present in fewer than N genes (default: 10)")
    parser.add_argument("--min-motifs-per-gene", type=int, default=2,
                        help="Genes with fewer motifs are excluded from clustering (default: 2)")
    parser.add_argument("--max-genes", type=int, default=5000,
                        help="Subsample down to this many genes for clustering (default: 5000)")
    parser.add_argument("--min-jaccard", type=float, default=0.1,
                        help="Edge kept if jaccard >= this (default: 0.1)")
    parser.add_argument("--min-lift", type=float, default=2.0,
                        help="Edge kept if lift >= this (default: 2.0)")
    parser.add_argument("--fdr", type=float, default=0.05,
                        help="BH q-value cutoff for edges and fingerprints (default: 0.05)")
    parser.add_argument("--n-clusters", type=int,
                        help="Force this many gene clusters (overrides --cluster-distance)")
    parser.add_argument("--cluster-distance", type=float, default=0.7,
                        help="Distance threshold for fcluster when --n-clusters absent (default: 0.7)")
    parser.add_argument("--skip-cooccurrence", action="store_true",
                        help="Skip motif co-occurrence analysis")
    parser.add_argument("--skip-clusters", action="store_true",
                        help="Skip gene clustering")
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args(argv)

    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)

    print("=== motif network / gene clustering ===")
    print(f"Loading {args.consolidated_file}...")
    df = load_consolidated(args.consolidated_file)
    print(f"  {len(df)} sites, {df['gene_id'].nunique()} genes, "
          f"{df['consolidated_motif_id'].nunique()} motifs")

    print(f"Building gene × motif presence matrix (>= {args.min_motif_sites} sites/motif)...")
    presence = build_presence_matrix(df, min_motif_sites=args.min_motif_sites)
    n_genes, n_motifs = presence.shape
    print(f"  matrix: {n_genes} genes × {n_motifs} motifs")
    if n_motifs < 2:
        print("ERROR: < 2 motifs survived filtering. Lower --min-motif-sites or "
              "provide more data.")
        return 1

    summary = [
        "# Motif network / gene clustering summary",
        f"Input: {args.consolidated_file}",
        f"Total sites: {len(df)} | genes: {df['gene_id'].nunique()} | "
        f"motifs: {df['consolidated_motif_id'].nunique()}",
        f"After filtering (>= {args.min_motif_sites} sites/motif): "
        f"{n_genes} genes × {n_motifs} motifs",
        "",
    ]

    if not args.skip_cooccurrence:
        print(f"Co-occurrence (jaccard>={args.min_jaccard} OR lift>={args.min_lift}; "
              f"q<={args.fdr})...")
        all_edges, kept = motif_cooccurrence(
            presence, args.min_jaccard, args.min_lift, args.fdr
        )
        edges_path = out / "motif_cooccurrence_edges.tsv"
        kept.to_csv(edges_path, sep="\t", index=False)
        print(f"  {len(kept)} / {len(all_edges)} edges passed -> {edges_path}")

        modules = connected_modules(kept, list(presence.columns))
        n_real_modules = len({v for v in modules.values() if v > 0})
        largest = max(
            (sum(1 for v in modules.values() if v == mid)
             for mid in {v for v in modules.values() if v > 0}),
            default=0,
        )
        modules_df = (pd.DataFrame([
            {"motif_id": m, "module_id": mid,
             "gene_count": int(presence[m].sum())}
            for m, mid in modules.items()
        ]).sort_values(["module_id", "gene_count"], ascending=[True, False]))
        modules_path = out / "motif_modules.tsv"
        modules_df.to_csv(modules_path, sep="\t", index=False)
        print(f"  {n_real_modules} multi-motif modules "
              f"(largest = {largest}) -> {modules_path}")

        gene_counts = presence.sum(axis=0).to_dict()
        nodes_meta = {m: (int(gene_counts.get(m, 0)), modules.get(m, 0))
                      for m in presence.columns}
        graphml_path = out / "motif_cooccurrence_network.graphml"
        write_graphml(kept, nodes_meta, graphml_path)
        print(f"  GraphML -> {graphml_path}")

        summary += [
            "## Motif co-occurrence",
            f"  Pairs tested: {len(all_edges)}",
            f"  Edges passing filters: {len(kept)}",
            f"  Multi-motif modules: {n_real_modules} (largest = {largest})",
            "",
        ]

    if not args.skip_clusters:
        print(f"Clustering genes (Jaccard avg-link; min {args.min_motifs_per_gene} "
              f"motifs/gene; max {args.max_genes} genes)...")
        cluster_df, sampled = cluster_genes(
            presence, min_motifs_per_gene=args.min_motifs_per_gene,
            max_genes=args.max_genes, n_clusters=args.n_clusters,
            dist_threshold=args.cluster_distance, random_state=args.seed,
        )
        clusters_path = out / "gene_clusters.tsv"
        cluster_df.to_csv(clusters_path, sep="\t", index=False)
        n_clustered = int(cluster_df["cluster_id"].notna().sum())
        n_actual = int(cluster_df["cluster_id"].dropna().nunique())
        print(f"  {n_clustered}/{n_genes} genes -> {n_actual} clusters"
              f"{' (subsampled)' if sampled else ''} -> {clusters_path}")

        print(f"Cluster motif fingerprints (Fisher, q<={args.fdr})...")
        fp = cluster_fingerprints(presence, cluster_df, fdr=args.fdr)
        fp_path = out / "gene_cluster_fingerprints.tsv"
        fp.to_csv(fp_path, sep="\t", index=False)
        print(f"  {len(fp)} cluster×motif enrichments -> {fp_path}")

        summary += [
            "## Gene clustering",
            f"  Clustered: {n_clustered} / {n_genes} genes in {n_actual} clusters",
            f"  Subsampled to {args.max_genes}? {sampled}",
            f"  Enriched cluster-motif pairs (q<={args.fdr}): {len(fp)}",
            "",
        ]

    pos = positional_summary(df, list(presence.columns))
    if not pos.empty:
        pos_path = out / "motif_positional_summary.tsv"
        pos.to_csv(pos_path, sep="\t", index=False)
        print(f"Positional summary -> {pos_path}")
        summary += ["## Positional summary",
                    f"  Motifs profiled: {len(pos)}", ""]

    summary_path = out / "network_summary.txt"
    summary_path.write_text("\n".join(summary))
    print(f"\nSummary -> {summary_path}")
    print(f"All outputs in {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
