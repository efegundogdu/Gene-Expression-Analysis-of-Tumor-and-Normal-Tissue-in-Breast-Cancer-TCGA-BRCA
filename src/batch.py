"""Stage 2: preprocessing, normalization, outlier flags and batch diagnostics."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.decomposition import PCA

SEED = 42

ROOT = Path(__file__).resolve().parents[1]
PROCESSED = ROOT / "data" / "processed"
TABLES = ROOT / "results" / "tables"

MIN_CPM = 1.0
MIN_PLATE_SIZE = 10
N_PCS = 10

# columns this stage adds to samples.csv
STAGE2_COLS = ["lib_size", "lib_size_pc", "size_factor", "batch_group", "batch_has_both",
               "lib_size_z", "PC1_z_within_condition", "PC2_z_within_condition", "outlier_flag"]


# ---------------------------------------------------------------- preprocessing

def load_inputs() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Return (counts restricted to kept samples, columns = aliquot barcode), genes, samples."""
    samples = pd.read_csv(PROCESSED / "samples.csv")
    samples = samples.drop(columns=[c for c in STAGE2_COLS if c in samples.columns])  # rerun-safe
    counts = pd.read_parquet(PROCESSED / "counts.parquet")
    genes = pd.read_parquet(PROCESSED / "genes.parquet")
    counts = counts[samples["file_id"]]
    counts.columns = samples["aliquot_barcode"].to_numpy()
    return counts, genes, samples.set_index("aliquot_barcode", drop=False)


def cpm(counts: pd.DataFrame) -> pd.DataFrame:
    return counts / counts.sum(axis=0) * 1e6


def filter_genes(counts: pd.DataFrame, genes: pd.DataFrame, n_min: int) -> pd.DataFrame:
    """Keep protein-coding genes with CPM >= 1 in at least n_min samples.

    CPM uses the protein-coding library size so the filter matches what is analysed.
    """
    pc = genes.index[genes["gene_type"] == "protein_coding"]
    sub = counts.loc[pc]
    keep = (cpm(sub) >= MIN_CPM).sum(axis=1) >= n_min
    return sub.loc[keep]


def size_factors(counts: pd.DataFrame) -> pd.Series:
    """DESeq2 median-of-ratios size factors (genes with a zero in any sample are skipped)."""
    with np.errstate(divide="ignore"):
        log_counts = np.log(counts.to_numpy(dtype=float))
    log_geo_mean = log_counts.mean(axis=1)
    usable = np.isfinite(log_geo_mean)
    ratios = log_counts[usable] - log_geo_mean[usable, None]
    sf = np.exp(np.median(ratios, axis=0))
    return pd.Series(sf, index=counts.columns, name="size_factor")


def log_cpm(counts: pd.DataFrame, prior: float = 1.0) -> pd.DataFrame:
    return np.log2(cpm(counts) + prior)


# ---------------------------------------------------------------- QC / PCA

def robust_z(x: pd.Series) -> pd.Series:
    mad = stats.median_abs_deviation(x, scale="normal")
    return (x - x.median()) / mad


def run_pca(logcpm: pd.DataFrame, n_pcs: int = N_PCS) -> tuple[pd.DataFrame, np.ndarray]:
    """PCA on samples using centred log-CPM of all filtered genes."""
    x = logcpm.T.to_numpy()
    x = x - x.mean(axis=0)
    pca = PCA(n_components=n_pcs, random_state=SEED, svd_solver="randomized")
    scores = pca.fit_transform(x)
    cols = [f"PC{i + 1}" for i in range(n_pcs)]
    return pd.DataFrame(scores, index=logcpm.columns, columns=cols), pca.explained_variance_ratio_


def flag_outliers(samples: pd.DataFrame, scores: pd.DataFrame,
                  z_lib: float = 3.0, z_pc: float = 4.0) -> pd.DataFrame:
    """Flag (not remove) samples with extreme library size or PC1/PC2 score.

    PC scores are judged within each condition, since PC1 is expected to separate
    tumor from normal.
    """
    out = pd.DataFrame(index=samples.index)
    out["lib_size_z"] = robust_z(np.log10(samples["lib_size"]))
    out["flag_lib_size"] = out["lib_size_z"].abs() > z_lib
    pc_flags = []
    for pc in ["PC1", "PC2"]:
        z = scores[pc].groupby(samples["condition"]).transform(robust_z)
        out[f"{pc}_z_within_condition"] = z
        pc_flags.append(z.abs() > z_pc)
    out["flag_pca"] = np.logical_or.reduce(pc_flags)
    out["outlier_flag"] = out["flag_lib_size"] | out["flag_pca"]
    return out


# ---------------------------------------------------------------- batch diagnostics

def r2_categorical(y: np.ndarray, groups: pd.Series) -> tuple[float, float]:
    """R2 and adjusted R2 of y regressed on a one-hot categorical (= eta squared)."""
    df = pd.DataFrame({"y": y, "g": groups.to_numpy()})
    fitted = df.groupby("g")["y"].transform("mean")
    ss_tot = ((df["y"] - df["y"].mean()) ** 2).sum()
    ss_res = ((df["y"] - fitted) ** 2).sum()
    r2 = 1 - ss_res / ss_tot
    n, k = len(df), df["g"].nunique() - 1
    adj = 1 - (1 - r2) * (n - 1) / (n - k - 1)
    return float(r2), float(adj)


def pc_variance_table(scores: pd.DataFrame, var_ratio: np.ndarray,
                      samples: pd.DataFrame, factors=("plate", "tss", "condition")) -> pd.DataFrame:
    rows = []
    for i, pc in enumerate(scores.columns):
        row = {"PC": pc, "var_explained": var_ratio[i]}
        for f in factors:
            r2, adj = r2_categorical(scores[pc].to_numpy(), samples.loc[scores.index, f])
            row[f"R2_{f}"] = r2
            row[f"adjR2_{f}"] = adj
        rows.append(row)
    return pd.DataFrame(rows)


def plate_condition_test(samples: pd.DataFrame, factor: str = "plate") -> dict:
    tab = pd.crosstab(samples[factor], samples["condition"])
    chi2, p, dof, expected = stats.chi2_contingency(tab)
    n = tab.to_numpy().sum()
    v = np.sqrt(chi2 / (n * (min(tab.shape) - 1)))
    return {"table": tab, "chi2": chi2, "p": p, "dof": dof, "cramers_v": v,
            "frac_expected_lt5": float((expected < 5).mean()),
            "tumor_only": sorted(tab.index[tab["normal"] == 0]),
            "normal_only": sorted(tab.index[tab["tumor"] == 0])}


# ---------------------------------------------------------------- plate grouping

def group_plates(samples: pd.DataFrame, min_size: int = MIN_PLATE_SIZE) -> pd.Series:
    """Merge small plates with their neighbours in processing order.

    TCGA plate IDs are issued sequentially (base-36: digits before letters), so the
    sorted plate ID is used as processing order. Plates are walked in order and
    accumulated until the running group has >= min_size samples; a small trailing
    group is merged into the previous one. Returns plate -> batch_group.
    """
    sizes = samples["plate"].value_counts().sort_index()
    groups, current, current_n = [], [], 0
    for plate, n in sizes.items():
        current.append(plate)
        current_n += n
        if current_n >= min_size:
            groups.append(current)
            current, current_n = [], 0
    if current:
        if groups:
            groups[-1].extend(current)
        else:
            groups.append(current)
    mapping = {}
    for i, members in enumerate(groups):
        name = f"B{i + 1:02d}_" + "+".join(members)
        for p in members:
            mapping[p] = name
    return pd.Series(mapping, name="batch_group")


# ---------------------------------------------------------------- main

def main() -> dict:
    TABLES.mkdir(parents=True, exist_ok=True)
    counts, genes, samples = load_inputs()
    n_normal = int((samples["condition"] == "normal").sum())

    samples["lib_size"] = counts.sum(axis=0)
    filt = filter_genes(counts, genes, n_min=n_normal)
    sf = size_factors(filt)
    norm = filt / sf
    samples["lib_size_pc"] = filt.sum(axis=0)
    samples["size_factor"] = sf

    filt.to_parquet(PROCESSED / "counts_filtered.parquet")
    norm.to_parquet(PROCESSED / "normcounts.parquet")

    lc = log_cpm(filt)
    scores, var_ratio = run_pca(lc)
    scores.to_parquet(PROCESSED / "pca_scores.parquet")
    pd.Series(var_ratio, index=scores.columns).to_csv(PROCESSED / "pca_var_ratio.csv")

    flags = flag_outliers(samples, scores)
    flags[flags["outlier_flag"]].join(samples[["condition", "plate", "lib_size"]]) \
        .to_csv(TABLES / "outlier_flags.csv")
    pcv = pc_variance_table(scores, var_ratio, samples)
    pcv.to_csv(TABLES / "pc_variance.csv", index=False)

    test_plate = plate_condition_test(samples, "plate")
    test_tss = plate_condition_test(samples, "tss")
    test_plate["table"].to_csv(TABLES / "plate_by_condition.csv")

    mapping = group_plates(samples)
    samples["batch_group"] = samples["plate"].map(mapping)
    bg = pd.crosstab(samples["batch_group"], samples["condition"])
    both = bg.index[(bg["tumor"] > 0) & (bg["normal"] > 0)]
    samples["batch_has_both"] = samples["batch_group"].isin(both)
    bg.assign(has_both=bg.index.isin(both)).to_csv(TABLES / "batch_groups.csv")

    samples = samples.join(flags[["lib_size_z", "PC1_z_within_condition",
                                  "PC2_z_within_condition", "outlier_flag"]])
    samples.reset_index(drop=True).to_csv(PROCESSED / "samples.csv", index=False)

    return {"counts_kept": counts, "filt": filt, "sf": sf, "samples": samples, "flags": flags,
            "pcv": pcv, "test_plate": test_plate, "test_tss": test_tss,
            "batch_table": bg, "both": both, "scores": scores, "var_ratio": var_ratio,
            "n_protein_coding": int((genes["gene_type"] == "protein_coding").sum())}


if __name__ == "__main__":
    res = main()
    s = res["samples"]
    print("protein-coding genes:", res["n_protein_coding"], "-> kept:", res["filt"].shape[0])
    print("size factors: min %.3f median %.3f max %.3f" % tuple(np.percentile(res["sf"], [0, 50, 100])))
    print("lib size (all genes, M): min %.1f median %.1f max %.1f" %
          tuple(np.percentile(s["lib_size"] / 1e6, [0, 50, 100])))
    fl = res["flags"]
    print("outliers: lib", int(fl["flag_lib_size"].sum()), "pca", int(fl["flag_pca"].sum()),
          "any", int(fl["outlier_flag"].sum()))
    print(s.loc[s["outlier_flag"], ["condition", "plate", "lib_size", "lib_size_z",
                                    "PC1_z_within_condition", "PC2_z_within_condition"]].to_string())
    print(res["pcv"].round(3).to_string(index=False))
    for k in ("test_plate", "test_tss"):
        t = res[k]
        print(k, "chi2=%.1f dof=%d p=%.3g V=%.3f frac_exp<5=%.2f" %
              (t["chi2"], t["dof"], t["p"], t["cramers_v"], t["frac_expected_lt5"]))
        print("  normal-only:", t["normal_only"], "| tumor-only n:", len(t["tumor_only"]))
    print("tumor-only plates:", res["test_plate"]["tumor_only"])
    print(res["batch_table"].to_string())
    print("batch groups:", len(res["batch_table"]), "with both conditions:", len(res["both"]))
    print("samples in both-condition batches:", int(s["batch_has_both"].sum()), "of", len(s),
          s[s["batch_has_both"]]["condition"].value_counts().to_dict())
