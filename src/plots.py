"""Figures. All saved as 300 dpi PNG under results/figures/."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap

ROOT = Path(__file__).resolve().parents[1]
FIGURES = ROOT / "results" / "figures"
DPI = 300

# Reference palette (light mode)
SERIES = ["#2a78d6", "#eb6834", "#1baf7a"]  # first three slots: safe for scatter
OTHER = "#c3c2b7"
INK = "#0b0b0b"
INK2 = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
SURFACE = "#fcfcfb"
SEQ_BLUE = ["#fcfcfb", "#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]
CONDITION_COLORS = {"normal": SERIES[0], "tumor": SERIES[1]}


def style() -> None:
    plt.rcParams.update({
        "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
        "axes.edgecolor": OTHER, "axes.labelcolor": INK2, "xtick.color": MUTED,
        "ytick.color": MUTED, "text.color": INK, "axes.grid": True, "grid.color": GRID,
        "grid.linewidth": 0.6, "axes.spines.top": False, "axes.spines.right": False,
        "font.family": ["Segoe UI", "DejaVu Sans"], "font.size": 9,
        "axes.titlesize": 10, "axes.titleweight": "bold", "legend.frameon": False,
    })


def save(fig, name: str) -> Path:
    FIGURES.mkdir(parents=True, exist_ok=True)
    path = FIGURES / name
    fig.savefig(path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    return path


def _scatter_highlight(ax, scores, labels: pd.Series, highlight: list[str], var_ratio, title,
                       pcs=("PC1", "PC2")):
    """Grey background for all samples; up to three highlighted categories on top."""
    px, py = pcs
    rest = ~labels.isin(highlight)
    if rest.any():
        ax.scatter(scores.loc[rest, px], scores.loc[rest, py], s=8, c=OTHER,
                   linewidths=0, label=f"other ({rest.sum()})", zorder=1)
    for color, cat in zip(SERIES, highlight):
        m = labels == cat
        ax.scatter(scores.loc[m, px], scores.loc[m, py], s=12, c=color,
                   edgecolors=SURFACE, linewidths=0.4, label=f"{cat} ({m.sum()})", zorder=2)
    ix, iy = int(px[2:]) - 1, int(py[2:]) - 1
    ax.set_xlabel(f"{px} ({var_ratio[ix]:.1%})")
    ax.set_ylabel(f"{py} ({var_ratio[iy]:.1%})")
    ax.set_title(title, loc="left")
    ax.legend(loc="best", fontsize=7, markerscale=1.2, handletextpad=0.2)


def fig2_pca(scores: pd.DataFrame, var_ratio, samples: pd.DataFrame, pcv: pd.DataFrame) -> Path:
    """PCA of log-CPM coloured by condition, plate and TSS (rows: PC1-PC2, PC1-PC3).

    Plates and TSS have ~40 levels, too many for distinguishable colours, so each
    panel highlights three levels and greys the rest; R2 per PC is in the title.
    """
    style()
    fig, axes = plt.subplots(2, 3, figsize=(13, 8.4))
    s = samples.loc[scores.index]
    top_plates = s["plate"].value_counts().index[:2].tolist()
    plates = ["A089"] + [p for p in top_plates if p != "A089"][:2]
    top_tss = s["tss"].value_counts().index[:3].tolist()
    # bottom row: plates with the lowest / highest mean PC3 (n >= 10) next to A089
    big = s.groupby("plate").filter(lambda g: len(g) >= 10)
    m3 = scores.loc[big.index, "PC3"].groupby(big["plate"]).mean().sort_values()
    plates_pc3 = ["A089", m3.index[0], m3.index[-1]]
    panels = [
        ("condition", ["normal", "tumor"], "Condition"),
        ("plate", plates, "Plate"),
        ("tss", top_tss, "Tissue source site"),
    ]
    letters = "ABCDEF"
    for row, pcs in enumerate([("PC1", "PC2"), ("PC1", "PC3")]):
        for col, (f, hl, name) in enumerate(panels):
            if row == 1 and f == "plate":
                hl = plates_pc3
            iy = int(pcs[1][2:]) - 1
            title = (f"{letters[row * 3 + col]}  {name}  (R² {pcs[0]} = {pcv.loc[0, f'R2_{f}']:.2f}, "
                     f"{pcs[1]} = {pcv.loc[iy, f'R2_{f}']:.2f})")
            _scatter_highlight(axes[row, col], scores, s[f], hl, var_ratio, title, pcs)
    fig.tight_layout()
    return save(fig, "fig2_pca.png")


def fig_pc_r2(pcv: pd.DataFrame) -> Path:
    """Supplementary: adjusted R2 of each PC on plate, TSS and condition."""
    style()
    fig, ax = plt.subplots(figsize=(7, 3.2))
    x = np.arange(len(pcv))
    w = 0.26
    for i, (f, lab) in enumerate([("condition", "condition"), ("plate", "plate"), ("tss", "TSS")]):
        ax.bar(x + (i - 1) * w, pcv[f"adjR2_{f}"].clip(lower=0), width=w - 0.03,
               color=SERIES[i], label=lab)
    ax.set_xticks(x, pcv["PC"])
    ax.set_ylabel("adjusted R²")
    ax.set_title("Variance of each PC explained by condition and batch", loc="left")
    ax.legend(ncols=3, loc="upper right")
    ax.grid(axis="x", visible=False)
    return save(fig, "figS_pc_r2.png")


def fig3_plate_heatmap(table: pd.DataFrame) -> Path:
    """Plate x condition sample counts, plates in processing order."""
    style()
    t = table[["normal", "tumor"]].sort_index()
    cmap = LinearSegmentedColormap.from_list("seqblue", SEQ_BLUE)
    fig, ax = plt.subplots(figsize=(3.4, 10))
    vals = t.to_numpy()
    im = ax.imshow(np.log1p(vals), aspect="auto", cmap=cmap)
    ax.grid(False)
    ax.set_xticks([0, 1], ["normal", "tumor"])
    ax.set_yticks(range(len(t)), t.index, fontsize=7)
    ax.tick_params(length=0)
    for sp in ax.spines.values():
        sp.set_visible(False)
    thresh = np.log1p(vals).max() * 0.55
    for i in range(vals.shape[0]):
        for j in range(vals.shape[1]):
            v = vals[i, j]
            ax.text(j, i, str(v), ha="center", va="center", fontsize=6.5,
                    color=(SURFACE if np.log1p(v) > thresh else INK2) if v else MUTED)
    ax.set_title("Samples per plate", loc="left")
    cb = fig.colorbar(im, ax=ax, fraction=0.08, pad=0.04)
    ticks = [0, 1, 5, 20, 90]
    cb.set_ticks(np.log1p(ticks), labels=[str(t) for t in ticks])
    cb.outline.set_visible(False)
    cb.ax.tick_params(colors=MUTED, labelsize=7)
    return save(fig, "fig3_plate_condition_heatmap.png")


def fig_library_size(samples: pd.DataFrame) -> Path:
    """Supplementary: library size per sample, flagged outliers marked."""
    style()
    fig, ax = plt.subplots(figsize=(7, 3))
    bins = np.linspace(samples["lib_size"].min(), samples["lib_size"].max(), 50) / 1e6
    for cond, z in [("tumor", 1), ("normal", 2)]:
        v = samples.loc[samples["condition"] == cond, "lib_size"] / 1e6
        ax.hist(v, bins=bins, density=True, histtype="step", linewidth=2,
                color=CONDITION_COLORS[cond], label=f"{cond} ({len(v)})", zorder=z)
    flagged = samples.loc[samples["outlier_flag"], "lib_size"] / 1e6
    ax.plot(flagged, np.zeros(len(flagged)), "|", color=INK, ms=10,
            label=f"flagged ({len(flagged)})", zorder=3)
    ax.set_xlabel("library size (million reads assigned to genes)")
    ax.set_ylabel("density")
    ax.set_title("Library size", loc="left")
    ax.legend()
    return save(fig, "figS_library_size.png")


def fig_pc3_by_plate(scores: pd.DataFrame, samples: pd.DataFrame) -> Path:
    """Supplementary: PC3 per plate in processing order (plate ID order)."""
    style()
    s = samples.loc[scores.index]
    order = sorted(s["plate"].unique())
    has_normal = s.groupby("plate")["condition"].apply(lambda x: (x == "normal").any())
    fig, ax = plt.subplots(figsize=(11, 3.4))
    data = [scores.loc[s["plate"] == p, "PC3"].to_numpy() for p in order]
    bp = ax.boxplot(data, positions=range(len(order)), widths=0.6, patch_artist=True,
                    showfliers=False, medianprops={"color": INK, "linewidth": 1})
    for patch, p in zip(bp["boxes"], order):
        patch.set_facecolor(SERIES[0] if has_normal[p] else OTHER)
        patch.set_edgecolor(SURFACE)
    for w in bp["whiskers"] + bp["caps"]:
        w.set_color(MUTED)
    ax.set_xticks(range(len(order)), order, rotation=90, fontsize=7)
    ax.set_ylabel("PC3 score")
    ax.set_title("PC3 by plate, in processing order (blue = plate contains normals)", loc="left")
    ax.grid(axis="x", visible=False)
    return save(fig, "figS_pc3_by_plate.png")


def make_stage2_figures(res: dict) -> list[Path]:
    s = res["samples"]
    return [fig2_pca(res["scores"], res["var_ratio"], s, res["pcv"]), fig_pc_r2(res["pcv"]),
            fig3_plate_heatmap(res["test_plate"]["table"]), fig_library_size(s),
            fig_pc3_by_plate(res["scores"], s)]


# ---------------------------------------------------------------- stage 3

def fig4_deg_counts(summary: pd.DataFrame) -> Path:
    """Horizontal bars: up and down DEG count per protocol."""
    style()
    s = summary.iloc[::-1]
    y = np.arange(len(s))
    h = 0.38
    fig, ax = plt.subplots(figsize=(7, 3.6))
    ax.barh(y + h / 2, s["n_up"], height=h - 0.04, color=CONDITION_COLORS["tumor"], label="up in tumor")
    ax.barh(y - h / 2, s["n_down"], height=h - 0.04, color=CONDITION_COLORS["normal"], label="down in tumor")
    for yi, (u, d) in enumerate(zip(s["n_up"], s["n_down"])):
        ax.text(u + 20, yi + h / 2, f"{u}", va="center", fontsize=7, color=INK2)
        ax.text(d + 20, yi - h / 2, f"{d}", va="center", fontsize=7, color=INK2)
    ax.set_yticks(y, s.index)
    ax.set_xlabel("DEGs (FDR < 0.05, |log2FC| > 1)")
    ax.set_title("DEG count per protocol", loc="left")
    ax.grid(axis="y", visible=False)
    ax.legend(loc="lower right")
    return save(fig, "fig4_deg_counts.png")


def fig_deg_jaccard(jac: pd.DataFrame) -> Path:
    """Supplementary: pairwise Jaccard similarity of DEG sets."""
    style()
    cmap = LinearSegmentedColormap.from_list("seqblue", SEQ_BLUE)
    fig, ax = plt.subplots(figsize=(4.6, 4))
    lo = float(np.floor(jac.to_numpy().min() * 10) / 10)
    im = ax.imshow(jac.to_numpy(), cmap=cmap, vmin=lo, vmax=1)
    ax.grid(False)
    ax.set_xticks(range(len(jac)), jac.columns)
    ax.set_yticks(range(len(jac)), jac.index)
    ax.tick_params(length=0)
    for sp in ax.spines.values():
        sp.set_visible(False)
    for i in range(len(jac)):
        for j in range(len(jac)):
            v = jac.iloc[i, j]
            ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=7,
                    color=SURFACE if (v - lo) / (1 - lo) > 0.45 else INK2)
    ax.set_title("DEG set Jaccard similarity", loc="left")
    cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cb.outline.set_visible(False)
    cb.ax.tick_params(colors=MUTED, labelsize=7)
    return save(fig, "figS_deg_jaccard.png")


# ---------------------------------------------------------------- stage 4

NULL_SERIES = [("FP_A", "Null A: shuffled labels"), ("FP_rand", "Null B-rand: tumors, random fake normals"),
               ("FP_conf", "Null B-conf: tumors, plate-confounded fake normals")]


def fig5_false_positives(summary: pd.DataFrame) -> Path:
    """FP per protocol under the three nulls, mean with 95% bootstrap CI; one panel per
    threshold. Symlog y axis so that near-zero rank-test values stay visible."""
    style()
    thresholds = list(dict.fromkeys(summary["threshold"]))
    fig, axes = plt.subplots(1, len(thresholds), figsize=(13, 4.4), sharey=True)
    for ax, thr in zip(np.atleast_1d(axes), thresholds):
        d = summary[summary["threshold"] == thr].set_index("protocol")
        x = np.arange(len(d))
        w = 0.27
        for i, (key, label) in enumerate(NULL_SERIES):
            mean = d[key].to_numpy()
            err = np.vstack([mean - d[f"{key}_ci_lo"].to_numpy(), d[f"{key}_ci_hi"].to_numpy() - mean])
            ax.bar(x + (i - 1) * w, mean, width=w - 0.03, color=SERIES[i], label=label, zorder=2)
            ax.errorbar(x + (i - 1) * w, mean, yerr=np.clip(err, 0, None), fmt="none",
                        ecolor=INK2, elinewidth=0.8, capsize=2, zorder=3)
        ax.set_yscale("symlog", linthresh=1)
        ax.set_yticks([0, 1, 10, 100, 1000], ["0", "1", "10", "100", "1000"])
        ax.set_ylim(0, 3000)
        ax.set_xticks(x, d.index)
        ax.grid(axis="x", visible=False)
        ax.set_title(f"{thr}", loc="left")
    np.atleast_1d(axes)[0].set_ylabel("false positives per null dataset (mean, 95% CI)")
    fig.tight_layout(rect=(0, 0.08, 1, 1))
    handles, labels = np.atleast_1d(axes)[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower left", ncols=3, fontsize=8, bbox_to_anchor=(0.01, 0.0))
    fig.text(0.99, 0.01, "P4, P6: full data (B = 100); P1, P2, P3, P5, P7: subsamples (B = 50)",
             fontsize=7.5, color=MUTED, ha="right")
    return save(fig, "fig5_false_positives.png")
