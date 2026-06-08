#!/usr/bin/env python3
"""
Create static visual summaries from motif network analysis outputs.

Expected inputs (from motif_network.py in a single directory):
  - motif_cooccurrence_edges.tsv
  - motif_modules.tsv
  - gene_clusters.tsv
  - gene_cluster_fingerprints.tsv (optional)

Outputs (PNG files in --output directory):
  - edge_strength_distributions.png
  - motif_degree_top.png
  - module_sizes_top.png
  - gene_cluster_sizes_top.png
  - cluster_fingerprint_heatmap.png (optional)
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def load_table(path):
    if path.exists():
        return pd.read_csv(path, sep="\t")
    return None


def ensure_output_dir(path):
    path.mkdir(parents=True, exist_ok=True)
    return path


def plot_edge_strength(edges, out_dir, dpi):
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5), constrained_layout=True)

    sns.histplot(edges["jaccard"], bins=30, kde=True, ax=axes[0], color="#1f77b4")
    axes[0].set_title("Jaccard Distribution")
    axes[0].set_xlabel("Jaccard")
    axes[0].set_ylabel("Edge count")

    sns.histplot(edges["lift"], bins=30, kde=True, ax=axes[1], color="#d62728")
    axes[1].set_title("Lift Distribution")
    axes[1].set_xlabel("Lift")
    axes[1].set_ylabel("Edge count")

    out = out_dir / "edge_strength_distributions.png"
    fig.savefig(out, dpi=dpi)
    plt.close(fig)
    return out


def plot_motif_degree(edges, out_dir, top_n, dpi):
    degree = pd.concat([
        edges[["motif_a"]].rename(columns={"motif_a": "motif_id"}),
        edges[["motif_b"]].rename(columns={"motif_b": "motif_id"}),
    ]).value_counts("motif_id").reset_index(name="degree")

    top = degree.head(top_n).sort_values("degree", ascending=True)

    fig, ax = plt.subplots(figsize=(8, max(4, 0.35 * len(top))), constrained_layout=True)
    sns.barplot(data=top, x="degree", y="motif_id", ax=ax, color="#2ca02c")
    ax.set_title(f"Top {len(top)} Motifs by Network Degree")
    ax.set_xlabel("Degree")
    ax.set_ylabel("Motif")

    out = out_dir / "motif_degree_top.png"
    fig.savefig(out, dpi=dpi)
    plt.close(fig)
    return out


def plot_module_sizes(modules, out_dir, top_n, dpi):
    mod = modules[modules["module_id"] > 0]
    sizes = mod.groupby("module_id", as_index=False).size().rename(columns={"size": "motif_count"})
    top = sizes.sort_values("motif_count", ascending=False).head(top_n).sort_values("motif_count", ascending=True)

    fig, ax = plt.subplots(figsize=(8, max(4, 0.35 * len(top))), constrained_layout=True)
    sns.barplot(data=top, x="motif_count", y="module_id", orient="h", ax=ax, color="#9467bd")
    ax.set_title(f"Top {len(top)} Motif Module Sizes")
    ax.set_xlabel("Motifs in module")
    ax.set_ylabel("Module ID")

    out = out_dir / "module_sizes_top.png"
    fig.savefig(out, dpi=dpi)
    plt.close(fig)
    return out


def plot_gene_cluster_sizes(gene_clusters, out_dir, top_n, dpi):
    valid = gene_clusters.dropna(subset=["cluster_id"]).copy()
    if valid.empty:
        return None

    valid["cluster_id"] = valid["cluster_id"].astype(int)
    sizes = valid.groupby("cluster_id", as_index=False).size().rename(columns={"size": "gene_count"})
    top = sizes.sort_values("gene_count", ascending=False).head(top_n).sort_values("gene_count", ascending=True)

    fig, ax = plt.subplots(figsize=(8, max(4, 0.35 * len(top))), constrained_layout=True)
    sns.barplot(data=top, x="gene_count", y="cluster_id", orient="h", ax=ax, color="#ff7f0e")
    ax.set_title(f"Top {len(top)} Gene Cluster Sizes")
    ax.set_xlabel("Genes in cluster")
    ax.set_ylabel("Cluster ID")

    out = out_dir / "gene_cluster_sizes_top.png"
    fig.savefig(out, dpi=dpi)
    plt.close(fig)
    return out


def plot_cluster_fingerprint_heatmap(fingerprints, out_dir, top_n, dpi):
    if fingerprints is None or fingerprints.empty:
        return None

    scored = fingerprints.copy()
    scored["score"] = -np.log10(scored["fisher_q"].clip(lower=1e-300))

    top_motifs = scored.groupby("motif_id")["score"].max().sort_values(ascending=False).head(top_n).index
    top_clusters = scored.groupby("cluster_id")["score"].max().sort_values(ascending=False).head(top_n).index

    sub = scored[scored["motif_id"].isin(top_motifs) & scored["cluster_id"].isin(top_clusters)]
    if sub.empty:
        return None

    mat = sub.pivot_table(index="cluster_id", columns="motif_id", values="score", fill_value=0.0)
    mat = mat.loc[sorted(mat.index)]

    fig, ax = plt.subplots(figsize=(max(8, 0.35 * mat.shape[1]), max(4, 0.4 * mat.shape[0])), constrained_layout=True)
    sns.heatmap(mat, cmap="mako", ax=ax)
    ax.set_title("Cluster-Motif Enrichment Heatmap (-log10 q)")
    ax.set_xlabel("Motif")
    ax.set_ylabel("Gene cluster")

    out = out_dir / "cluster_fingerprint_heatmap.png"
    fig.savefig(out, dpi=dpi)
    plt.close(fig)
    return out


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Create static figures from motif network outputs."
    )
    parser.add_argument("network_dir", help="Directory with motif_network.py outputs")
    parser.add_argument("--output", "-o", help="Output directory for figures (default: <network_dir>/figures)")
    parser.add_argument("--top-n", type=int, default=20, help="Top N categories to show in bar/heatmap plots")
    parser.add_argument("--dpi", type=int, default=200, help="Figure DPI (default: 200)")
    args = parser.parse_args(argv)

    network_dir = Path(args.network_dir)
    if not network_dir.exists():
        print(f"Error: network directory not found: {network_dir}")
        return 1

    out_dir = Path(args.output) if args.output else network_dir / "figures"
    ensure_output_dir(out_dir)

    sns.set_theme(style="whitegrid", context="talk")

    edges = load_table(network_dir / "motif_cooccurrence_edges.tsv")
    modules = load_table(network_dir / "motif_modules.tsv")
    gene_clusters = load_table(network_dir / "gene_clusters.tsv")
    fingerprints = load_table(network_dir / "gene_cluster_fingerprints.tsv")

    produced = []

    if edges is not None and not edges.empty and {"motif_a", "motif_b", "jaccard", "lift"}.issubset(edges.columns):
        produced.append(plot_edge_strength(edges, out_dir, args.dpi))
        produced.append(plot_motif_degree(edges, out_dir, args.top_n, args.dpi))
    else:
        print("[warn] Skipping edge plots: motif_cooccurrence_edges.tsv missing or malformed")

    if modules is not None and not modules.empty and {"module_id", "motif_id"}.issubset(modules.columns):
        produced.append(plot_module_sizes(modules, out_dir, args.top_n, args.dpi))
    else:
        print("[warn] Skipping module size plot: motif_modules.tsv missing or malformed")

    if gene_clusters is not None and not gene_clusters.empty and {"cluster_id", "gene_id"}.issubset(gene_clusters.columns):
        out = plot_gene_cluster_sizes(gene_clusters, out_dir, args.top_n, args.dpi)
        if out is not None:
            produced.append(out)
    else:
        print("[warn] Skipping gene cluster plot: gene_clusters.tsv missing or malformed")

    if fingerprints is not None and not fingerprints.empty and {"cluster_id", "motif_id", "fisher_q"}.issubset(fingerprints.columns):
        out = plot_cluster_fingerprint_heatmap(fingerprints, out_dir, args.top_n, args.dpi)
        if out is not None:
            produced.append(out)
    else:
        print("[warn] Skipping fingerprint heatmap: gene_cluster_fingerprints.tsv missing or malformed")

    produced = [p for p in produced if p is not None]
    if not produced:
        print("No figures were created. Check that network outputs are present.")
        return 1

    print("Created figures:")
    for fig in produced:
        print(f"  - {fig}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
