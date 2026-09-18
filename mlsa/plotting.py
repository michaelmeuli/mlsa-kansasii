"""Visualizations shared by main1_locus_discovery and
main2_sanger_differentiation: an alignment-difference heatmap (shows exactly
where sequences differ) and a species-colored tree image. Both are saved as
PNG+SVG using one consistent species-color mapping (mlsa.SPECIES_COLORS)
across the whole repo.
"""
from __future__ import annotations

import collections
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch

from . import SPECIES_COLORS

_MATCH, _MISMATCH, _GAP = 0, 1, 2
_HEATMAP_CMAP = ListedColormap(["#f7f7f7", "#b2182b", "#999999"])


def plot_alignment_heatmap(alignment: dict[str, str], species_of: dict[str, str],
                            out_prefix: Path, title: str) -> None:
    """Rows = sequences (sorted/colored by species), columns = variable
    alignment positions only (constant columns carry no information about
    differences and would just waste width), colored by whether each base
    matches the per-column majority consensus, mismatches it, or is a gap."""
    names = list(alignment.keys())
    seqs = [alignment[n] for n in names]
    length = len(seqs[0]) if seqs else 0

    var_cols = [
        col for col in range(length)
        if len({s[col].upper() for s in seqs if s[col] != "-"}) > 1
    ]
    if not var_cols:
        var_cols = list(range(length))

    consensus = []
    for col in var_cols:
        counts = collections.Counter(s[col].upper() for s in seqs if s[col] != "-")
        consensus.append(counts.most_common(1)[0][0] if counts else "N")

    mat = np.full((len(names), len(var_cols)), _MATCH, dtype=int)
    for i, s in enumerate(seqs):
        for j, col in enumerate(var_cols):
            c = s[col].upper()
            if c == "-":
                mat[i, j] = _GAP
            elif c != consensus[j]:
                mat[i, j] = _MISMATCH

    order = sorted(range(len(names)), key=lambda i: (species_of.get(names[i], "~"), names[i]))
    names = [names[i] for i in order]
    mat = mat[order]

    fig_h = max(3.0, 0.18 * len(names) + 1.5)
    fig_w = min(max(6.0, 0.04 * len(var_cols) + 3.0), 22.0)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    ax.imshow(mat, aspect="auto", cmap=_HEATMAP_CMAP, vmin=0, vmax=2, interpolation="nearest")
    ax.set_yticks(range(len(names)))
    ax.set_yticklabels(names, fontsize=6)
    for tick_label, name in zip(ax.get_yticklabels(), names):
        sp = species_of.get(name)
        if sp in SPECIES_COLORS:
            tick_label.set_color(SPECIES_COLORS[sp])
    ax.set_xticks([])
    ax.set_xlabel(f"variable alignment positions (n={len(var_cols)} of {length} total)")
    ax.set_title(title)

    state_handles = [
        Patch(facecolor="#f7f7f7", edgecolor="black", label="matches consensus"),
        Patch(facecolor="#b2182b", label="differs from consensus"),
        Patch(facecolor="#999999", label="gap"),
    ]
    ax.legend(handles=state_handles, loc="upper center", bbox_to_anchor=(0.5, -0.12),
              ncol=3, frameon=False, fontsize=7)

    present_species = {sp for sp in species_of.values() if sp in SPECIES_COLORS}
    species_handles = [Patch(facecolor=SPECIES_COLORS[sp], label=sp) for sp in sorted(present_species)]
    if species_handles:
        fig.legend(handles=species_handles, loc="upper left", bbox_to_anchor=(1.0, 1.0),
                   fontsize=7, title="species", frameon=False)

    fig.tight_layout()
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(f"{out_prefix}.png", dpi=200, bbox_inches="tight")
    fig.savefig(f"{out_prefix}.svg", bbox_inches="tight")
    plt.close(fig)


def plot_tree(treefile: Path, species_of: dict[str, str], out_prefix: Path, title: str) -> None:
    """Render a Newick tree with species-colored tip labels using Biopython's
    matplotlib-backed Phylo.draw (avoids the heavier ete3/Qt dependency)."""
    from Bio import Phylo

    tree = Phylo.read(str(treefile), "newick")
    n_tips = tree.count_terminals()
    fig_h = max(3.0, 0.28 * n_tips + 1.5)
    fig, ax = plt.subplots(figsize=(8.0, fig_h))

    def label_colors(name: str) -> str:
        return SPECIES_COLORS.get(species_of.get(name, ""), "black")

    Phylo.draw(tree, axes=ax, do_show=False, label_colors=label_colors)
    ax.set_title(title)
    fig.tight_layout()
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(f"{out_prefix}.png", dpi=200, bbox_inches="tight")
    fig.savefig(f"{out_prefix}.svg", bbox_inches="tight")
    plt.close(fig)
