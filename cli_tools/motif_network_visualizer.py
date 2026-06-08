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
from matplotlib.collections import LineCollection
from matplotlib import cm


def load_table(path):
    if path.exists():
        return pd.read_csv(path, sep="\t")
    return None


def ensure_output_dir(path):
    path.mkdir(parents=True, exist_ok=True)
    return path


def compute_degree(edges):
    return pd.concat([
        edges[["motif_a"]].rename(columns={"motif_a": "motif_id"}),
        edges[["motif_b"]].rename(columns={"motif_b": "motif_id"}),
    ]).value_counts("motif_id").reset_index(name="degree")


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
    degree = compute_degree(edges)

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


def _plot_network(edges, modules, out_file, title, max_nodes, label_top, dpi):
    degree = compute_degree(edges)
    nodes = degree.head(max_nodes)["motif_id"].tolist()
    if not nodes:
        return None

    sub = edges[edges["motif_a"].isin(nodes) & edges["motif_b"].isin(nodes)].copy()
    if sub.empty:
        return None

    degree_map = dict(zip(degree["motif_id"], degree["degree"]))
    module_map = {}
    if modules is not None and not modules.empty and {"motif_id", "module_id"}.issubset(modules.columns):
        module_map = dict(zip(modules["motif_id"], modules["module_id"]))

    node_df = pd.DataFrame({"motif_id": nodes})
    node_df["degree"] = node_df["motif_id"].map(degree_map).fillna(0)
    node_df["module_id"] = node_df["motif_id"].map(module_map).fillna(0).astype(int)
    node_df = node_df.sort_values(["module_id", "degree"], ascending=[True, False]).reset_index(drop=True)

    n = len(node_df)
    angles = np.linspace(0, 2 * np.pi, n, endpoint=False)
    node_df["x"] = np.cos(angles)
    node_df["y"] = np.sin(angles)
    pos = dict(zip(node_df["motif_id"], zip(node_df["x"], node_df["y"])))

    segments = []
    strengths = []
    for _, row in sub.iterrows():
        a = row["motif_a"]
        b = row["motif_b"]
        if a not in pos or b not in pos:
            continue
        segments.append([pos[a], pos[b]])
        strengths.append(float(row.get("jaccard", 0.0)))

    fig, ax = plt.subplots(figsize=(10, 10), constrained_layout=True)
    ax.set_title(title)

    if segments:
        strengths_arr = np.array(strengths)
        if np.allclose(strengths_arr.max(), strengths_arr.min()):
            widths = np.full_like(strengths_arr, 0.8)
        else:
            widths = 0.4 + 2.6 * (strengths_arr - strengths_arr.min()) / (strengths_arr.max() - strengths_arr.min())
        edge_collection = LineCollection(segments, colors="#7f8c8d", linewidths=widths, alpha=0.25, zorder=1)
        ax.add_collection(edge_collection)

    modules_unique = sorted(node_df["module_id"].unique())
    color_lookup = {}
    palette = cm.get_cmap("tab20", max(1, len([m for m in modules_unique if m > 0])))
    idx = 0
    for mid in modules_unique:
        if mid <= 0:
            color_lookup[mid] = "#bdc3c7"
        else:
            color_lookup[mid] = palette(idx)
            idx += 1

    node_sizes = 30 + 8 * np.sqrt(node_df["degree"].values)
    node_colors = [color_lookup[m] for m in node_df["module_id"]]
    ax.scatter(node_df["x"], node_df["y"], s=node_sizes, c=node_colors,
               edgecolor="white", linewidth=0.5, zorder=3)

    if label_top > 0:
        top_labels = node_df.nlargest(label_top, "degree")
        for _, row in top_labels.iterrows():
            ax.text(row["x"] * 1.08, row["y"] * 1.08, row["motif_id"],
                    fontsize=8, ha="center", va="center", zorder=4)

    ax.set_aspect("equal")
    ax.set_xlim(-1.25, 1.25)
    ax.set_ylim(-1.25, 1.25)
    ax.axis("off")

    fig.savefig(out_file, dpi=dpi)
    plt.close(fig)
    return out_file


def plot_network_overview(edges, modules, out_dir, max_nodes, label_top, dpi):
    return _plot_network(
        edges,
        modules,
        out_dir / "motif_network_overview.png",
        f"Motif Co-occurrence Network (Top {max_nodes} nodes)",
        max_nodes,
        label_top,
        dpi,
    )


def plot_network_main_module(edges, modules, out_dir, max_nodes, label_top, dpi):
    if modules is None or modules.empty or "module_id" not in modules.columns:
        return None

    non_singleton = modules[modules["module_id"] > 0]
    if non_singleton.empty:
        return None
    main_module = int(non_singleton["module_id"].value_counts().idxmax())
    main_nodes = set(non_singleton.loc[non_singleton["module_id"] == main_module, "motif_id"])

    sub_edges = edges[edges["motif_a"].isin(main_nodes) & edges["motif_b"].isin(main_nodes)].copy()
    if sub_edges.empty:
        return None

    return _plot_network(
        sub_edges,
        modules,
        out_dir / "motif_network_main_module.png",
        f"Motif Co-occurrence Network (Main module {main_module})",
        max_nodes,
        label_top,
        dpi,
    )


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
    parser.add_argument("--network-nodes", type=int, default=120,
                        help="Max nodes shown in network graphs (default: 120)")
    parser.add_argument("--label-top", type=int, default=15,
                        help="Label top N high-degree nodes in network graphs (default: 15)")
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
        produced.append(plot_network_overview(
            edges, modules, out_dir, args.network_nodes, args.label_top, args.dpi
        ))
        produced.append(plot_network_main_module(
            edges, modules, out_dir, args.network_nodes, args.label_top, args.dpi
        ))
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
