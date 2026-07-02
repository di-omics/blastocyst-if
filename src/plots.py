#!/usr/bin/env python3
"""QC plots for the blastocyst IF pipeline (writes assets/blastocyst_qc.png).

Three panels: nucleus count + volume distribution, marker intensity
distributions by lineage, and classification proportions.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import omics_style as S; S.apply()

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUTPUTS = ROOT / "outputs"
ASSETS = ROOT / "assets"; ASSETS.mkdir(exist_ok=True)

LINEAGE_COLORS = {"EPI": "blue", "PE": "green", "TE": "peach"}


def main():
    classified = pd.read_csv(OUTPUTS / "classified.tsv", sep="\t")
    props = pd.read_csv(OUTPUTS / "proportions.tsv", sep="\t")

    fig, ax = plt.subplots(1, 3, figsize=(12, 3.8))

    # (0) nucleus volume distribution by lineage
    lineages = props.lineage.tolist()
    for lin in lineages:
        sub = classified[classified.lineage == lin]
        key = LINEAGE_COLORS.get(lin, "lav")
        ax[0].hist(sub.volume_voxels, bins=15, alpha=0.55, label=lin,
                   color=S.PALETTE[key], edgecolor=S.OUTLINE[key], linewidth=0.8)
    ax[0].set(title=f"Nucleus volume ({len(classified)} nuclei)",
              xlabel="volume (voxels)", ylabel="count")
    ax[0].legend()

    # (1) marker intensity distributions by lineage
    markers = [c for c in classified.columns if c.endswith("_mean") and c != "DAPI_mean"]
    x_pos = np.arange(len(markers))
    bar_w = 0.8 / len(lineages)
    for i, lin in enumerate(lineages):
        sub = classified[classified.lineage == lin]
        means = [sub[m].mean() for m in markers]
        key = LINEAGE_COLORS.get(lin, "lav")
        ax[1].bar(x_pos + i * bar_w, means, width=bar_w, label=lin,
                  color=S.PALETTE[key], edgecolor=S.OUTLINE[key], linewidth=0.8)
    ax[1].set(title="Mean marker intensity by lineage",
              ylabel="mean intensity", xticks=x_pos + bar_w)
    ax[1].set_xticklabels([m.replace("_mean", "") for m in markers], fontsize=8)
    ax[1].legend()
    ax[1].grid(axis="x", visible=False)

    # (2) classification proportions
    colors = [S.PALETTE[LINEAGE_COLORS.get(ln, "lav")] for ln in lineages]
    outlines = [S.OUTLINE[LINEAGE_COLORS.get(ln, "lav")] for ln in lineages]
    ax[2].bar(lineages, props.fraction * 100, color=colors,
              edgecolor=outlines, linewidth=0.9)
    ax[2].set(title="Lineage proportions", ylabel="% of nuclei")
    ax[2].grid(axis="x", visible=False)
    for i, row in props.iterrows():
        ax[2].text(i, row.fraction * 100 + 1.5, f"{row['count']}",
                   ha="center", fontsize=9, color=S.MUTED)

    out = ASSETS / "blastocyst_qc.png"
    fig.savefig(out)
    print("wrote", out)


if __name__ == "__main__":
    main()
