"""Stage 5: benchmark of the DE protocols against a batch-free paired reference.

Reference: patients whose tumor and normal sit on the same plate (no plate confounding
within a pair). (a) Wilcoxon signed-rank on within-patient differences of
log2(normalized count + 1); (b) PyDESeq2 ~ patient + condition. A gene is in the
reference if it is a DEG (FDR < 0.05, |log2FC| > 1) in both, with the same direction.
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
import de_methods as dm  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
PROCESSED = ROOT / "data" / "processed"
TABLES = ROOT / "results" / "tables"
MIN_PAIRS = 10


def paired_samples(samples: pd.DataFrame) -> pd.DataFrame:
    """same_plate_paired samples, sorted so tumor/normal rows line up per patient."""
    s = samples[samples["same_plate_paired"]].copy()
    return s.sort_values(["patient", "condition"])


def reference_wilcoxon(counts: pd.DataFrame, s: pd.DataFrame) -> pd.DataFrame:
    norm = dm.normalize(counts[s.index])
    lg = np.log2(norm + 1)
    t = lg[s.index[s["condition"] == "tumor"]].to_numpy()
    n = lg[s.index[s["condition"] == "normal"]].to_numpy()
    # columns line up by patient (both sorted by patient)
    pt_t = s.loc[s["condition"] == "tumor", "patient"].to_numpy()
    pt_n = s.loc[s["condition"] == "normal", "patient"].to_numpy()
    assert (pt_t == pt_n).all()
    d = t - n
    res = stats.wilcoxon(d, axis=1, alternative="two-sided", zero_method="wilcox",
                         method="approx")
    p = np.asarray(res.pvalue, dtype=float)
    out = pd.DataFrame({"log2FC": np.median(d, axis=1), "pvalue": p}, index=counts.index)
    return dm.finalize(out)


def reference_deseq(counts: pd.DataFrame, s: pd.DataFrame, n_cpus: int | None = None) -> pd.DataFrame:
    meta = pd.DataFrame({"condition": pd.Categorical(s["condition"], categories=["normal", "tumor"]),
                         "patient": s["patient"].to_numpy()}, index=s.index)
    return dm.run_pydeseq2(counts[s.index], meta, "~patient + condition", n_cpus)


def build_reference(n_cpus: int | None = None) -> dict:
    counts, samples, genes = dm.load_data()
    s = paired_samples(samples)
    n_pairs = int(s["patient"].nunique())
    if n_pairs < MIN_PAIRS:
        raise RuntimeError(f"only {n_pairs} same-plate pairs (< {MIN_PAIRS}); ask before continuing")
    t0 = time.perf_counter()
    w = reference_wilcoxon(counts, s)
    t_w = time.perf_counter() - t0
    t0 = time.perf_counter()
    d = reference_deseq(counts, s, n_cpus)
    t_d = time.perf_counter() - t0
    up = (w["direction"] == "up") & (d["direction"] == "up")
    down = (w["direction"] == "down") & (d["direction"] == "down")
    ref = pd.DataFrame({"gene_name": genes.loc[counts.index, "gene_name"].to_numpy(),
                        "wilcoxon_log2FC": w["log2FC"], "wilcoxon_FDR": w["FDR"],
                        "wilcoxon_DEG": w["direction"],
                        "deseq_log2FC": d["log2FC"], "deseq_FDR": d["FDR"],
                        "deseq_DEG": d["direction"]}, index=counts.index)
    ref["reference"] = np.where(up, "up", np.where(down, "down", ""))
    ref.index.name = "gene_id"
    ref.to_csv(TABLES / "reference_same_plate_paired.csv")
    info = {"n_pairs": n_pairs, "n_samples": len(s), "n_plates": int(s["plate"].nunique()),
            "wilcoxon_DEG": int(w["is_DEG"].sum()), "deseq_DEG": int(d["is_DEG"].sum()),
            "both_any_direction": int((w["is_DEG"] & d["is_DEG"]).sum()),
            "direction_conflicts": int((w["is_DEG"] & d["is_DEG"] & (w["direction"] != d["direction"])).sum()),
            "reference_up": int(up.sum()), "reference_down": int(down.sum()),
            "reference_total": int((up | down).sum()),
            "seconds_wilcoxon": t_w, "seconds_deseq": t_d}
    (TABLES / "reference_info.json").write_text(json.dumps(info, indent=1))
    return {"reference": ref, "info": info}


def load_reference() -> pd.DataFrame:
    ref = pd.read_csv(TABLES / "reference_same_plate_paired.csv", index_col="gene_id")
    for c in ["reference", "wilcoxon_DEG", "deseq_DEG"]:
        ref[c] = ref[c].fillna("")  # empty strings are read back as NaN
    return ref


def protocol_deg(p: str) -> pd.DataFrame:
    return dm.load_protocol(p)


def benchmark(ref: pd.DataFrame, protocols=dm.PROTOCOLS) -> pd.DataFrame:
    """Precision / recall / Jaccard of each protocol's DEG list vs the reference.

    A protocol DEG counts as a true positive only if it is in the reference with the
    same direction.
    """
    ref_set = {(g, d) for g, d in ref["reference"].items() if d in ("up", "down")}
    rows = []
    for p in protocols:
        r = protocol_deg(p)
        deg = {(g, d) for g, d in r.loc[r["is_DEG"], "direction"].items()}
        tp = len(deg & ref_set)
        rows.append({"protocol": p, "n_DEG": len(deg), "n_reference": len(ref_set), "TP": tp,
                     "precision": tp / len(deg) if deg else np.nan,
                     "recall": tp / len(ref_set) if ref_set else np.nan,
                     "jaccard": tp / len(deg | ref_set) if (deg | ref_set) else np.nan,
                     "wrong_direction": len({g for g, _ in deg} & {g for g, _ in ref_set}) - tp})
    out = pd.DataFrame(rows).set_index("protocol")
    out.to_csv(TABLES / "benchmark.csv")
    return out


def marker_table(ref: pd.DataFrame, protocols=dm.PROTOCOLS) -> pd.DataFrame:
    rows = []
    for g, expected in dm.MARKERS.items():
        row = {"gene": g, "expected": expected}
        rr = ref[ref["gene_name"] == g]
        row["reference"] = rr["reference"].iloc[0] if len(rr) else "not tested"
        row["reference_log2FC"] = rr["deseq_log2FC"].iloc[0] if len(rr) else np.nan
        for p in protocols:
            r = protocol_deg(p)
            h = r[r["gene_name"] == g]
            row[f"{p}_log2FC"] = h["log2FC"].iloc[0]
            row[f"{p}_FDR"] = h["FDR"].iloc[0]
            row[f"{p}_found"] = bool(h["is_DEG"].iloc[0])
            row[f"{p}_direction_ok"] = bool(h["is_DEG"].iloc[0]) and h["direction"].iloc[0] == expected
        rows.append(row)
    out = pd.DataFrame(rows).set_index("gene")
    out.to_csv(TABLES / "benchmark_markers.csv")
    return out


def deg_sets(protocols=dm.PROTOCOLS) -> dict[str, set]:
    sets = {}
    for p in protocols:
        r = protocol_deg(p)
        sets[p] = set(r.index[r["is_DEG"]])
    return sets


def benchmark_common_lfc(ref: pd.DataFrame, protocols=dm.PROTOCOLS) -> pd.DataFrame:
    """Sensitivity: the same benchmark with the common-log2FC DEG lists (Stage 3)."""
    ref_set = {(g, d) for g, d in ref["reference"].items() if d in ("up", "down")}
    rows = []
    for p in protocols:
        c = pd.read_csv(TABLES / "de_common_lfc" / f"de_{p}_common_lfc.csv", index_col=0)
        c["direction"] = c["direction"].fillna("")
        deg = {(g, d) for g, d in c.loc[c["is_DEG"], "direction"].items()}
        tp = len(deg & ref_set)
        rows.append({"protocol": p, "n_DEG": len(deg), "TP": tp, "precision": tp / len(deg),
                     "recall": tp / len(ref_set), "jaccard": tp / len(deg | ref_set)})
    out = pd.DataFrame(rows).set_index("protocol")
    out.to_csv(TABLES / "benchmark_common_lfc.csv")
    return out


def run_all(rebuild_reference: bool = False) -> dict:
    import plots
    if rebuild_reference or not (TABLES / "reference_same_plate_paired.csv").exists():
        build_reference()
    ref = load_reference()
    res = {"reference_info": json.loads((TABLES / "reference_info.json").read_text()),
           "benchmark": benchmark(ref), "benchmark_common_lfc": benchmark_common_lfc(ref),
           "markers": marker_table(ref)}
    res["fig7"] = plots.fig7_markers(res["markers"], dm.PROTOCOLS)
    res["fig6"] = plots.fig6_upset(deg_sets())
    return res


if __name__ == "__main__":
    import warnings
    warnings.simplefilter("ignore")
    out = run_all()
    print(json.dumps(out["reference_info"], indent=1))
    print(out["benchmark"].round(3).to_string())
    print(out["benchmark_common_lfc"].round(3).to_string())
