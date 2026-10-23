"""Stage 7: Table 1 (results/tables/summary.csv), one row per protocol."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import de_methods as dm  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
TABLES = ROOT / "results" / "tables"

DESCRIPTION = {
    "P1": "PyDESeq2 ~condition (apeGLM-shrunken log2FC)",
    "P2": "PyDESeq2 ~batch_group + condition",
    "P3": "ComBat-seq, then PyDESeq2 ~condition",
    "P4": "Wilcoxon rank-sum on normalized counts",
    "P5": "ComBat-seq, then Wilcoxon rank-sum",
    "P6": "van Elteren stratified Wilcoxon (strata = batch_group)",
    "P7": "Winsorized counts (95th pct, within condition), PyDESeq2 ~condition",
}


def build_table1() -> pd.DataFrame:
    info = pd.read_csv(TABLES / "de_protocol_summary.csv", index_col="protocol")
    perm = pd.read_csv(TABLES / "perm_summary.csv")
    bench = pd.read_csv(TABLES / "benchmark.csv", index_col="protocol")
    markers = pd.read_csv(TABLES / "benchmark_markers.csv", index_col="gene")
    rows = []
    for p in dm.PROTOCOLS:
        a = perm[(perm["protocol"] == p) & (perm["threshold"] == "FDR<0.05")].iloc[0]
        b = perm[(perm["protocol"] == p) & (perm["threshold"] != "FDR<0.05")].iloc[0]
        rows.append({
            "protocol": p, "method": DESCRIPTION[p],
            "n_samples_informative": int(info.loc[p, "n_used"]),
            "n_DEG": int(info.loc[p, "n_DEG"]), "n_up": int(info.loc[p, "n_up"]),
            "n_down": int(info.loc[p, "n_down"]),
            # false positives per null dataset (mean over replicates), FDR < 0.05
            "FP_A": a["FP_A"], "FP_rand": a["FP_rand"], "FP_conf": a["FP_conf"],
            "FP_batch": a["FP_batch"], "FP_batch_ci_lo": a["FP_batch_ci_lo"],
            "FP_batch_ci_hi": a["FP_batch_ci_hi"], "FP_C": a["FP_C"], "FDR_A": a["FDR_A"],
            # same with FDR < 0.05 and |log2FC| > 1
            "FP_A_lfc": b["FP_A"], "FP_rand_lfc": b["FP_rand"], "FP_conf_lfc": b["FP_conf"],
            "FP_batch_lfc": b["FP_batch"], "FDR_A_lfc": b["FDR_A"],
            "null_B_replicates": int(a["FP_conf_B"]), "null_dataset_A": a["FP_A_dataset"],
            "precision_vs_reference": bench.loc[p, "precision"],
            "recall_vs_reference": bench.loc[p, "recall"],
            "jaccard_vs_reference": bench.loc[p, "jaccard"],
            "markers_found_5": int(markers[f"{p}_direction_ok"].sum()),
            "runtime_s_full_data": info.loc[p, "seconds"],
        })
    t1 = pd.DataFrame(rows).set_index("protocol")
    t1.to_csv(TABLES / "summary.csv")
    return t1


if __name__ == "__main__":
    pd.set_option("display.width", 250)
    t = build_table1()
    print(t.drop(columns=["method"]).round(3).T.to_string())
