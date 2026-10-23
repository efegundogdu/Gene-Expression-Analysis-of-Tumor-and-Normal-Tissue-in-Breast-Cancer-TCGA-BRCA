"""Stage 3: differential expression protocols P1-P7.

Every protocol returns a table indexed by gene_id with:
log2FC, pvalue, FDR (BH over all tested genes), is_DEG, direction.
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import tempfile
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

SEED = 42
FDR_CUT = 0.05
LFC_CUT = 1.0

ROOT = Path(__file__).resolve().parents[1]
PROCESSED = ROOT / "data" / "processed"
TABLES = ROOT / "results" / "tables"
PROTOCOLS = ["P1", "P2", "P3", "P4", "P5", "P6", "P7"]
DESEQ_PROTOCOLS = {"P1", "P2", "P3", "P7"}
COMBAT_FILE = PROCESSED / "counts_combatseq.parquet"
COEF = "condition[T.tumor]"


# ---------------------------------------------------------------- shared helpers

def bh(p: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg over non-missing p-values; NaN stays NaN."""
    p = np.asarray(p, dtype=float)
    out = np.full_like(p, np.nan)
    ok = ~np.isnan(p)
    if ok.any():
        out[ok] = multipletests(p[ok], method="fdr_bh")[1]
    return out


def finalize(df: pd.DataFrame) -> pd.DataFrame:
    """Add FDR, DEG call and direction from log2FC and pvalue columns."""
    df = df.copy()
    df["FDR"] = bh(df["pvalue"].to_numpy())
    df["is_DEG"] = (df["FDR"] < FDR_CUT) & (df["log2FC"].abs() > LFC_CUT)
    df["direction"] = np.where(df["is_DEG"], np.where(df["log2FC"] > 0, "up", "down"), "")
    return df


def size_factors(counts: np.ndarray) -> np.ndarray:
    """DESeq2 median-of-ratios; counts is genes x samples."""
    with np.errstate(divide="ignore"):
        lc = np.log(counts.astype(float))
    lg = lc.mean(axis=1)
    ok = np.isfinite(lg)
    return np.exp(np.median(lc[ok] - lg[ok, None], axis=0))


def normalize(counts: pd.DataFrame) -> pd.DataFrame:
    return counts / size_factors(counts.to_numpy())


def median_log2fc(norm: np.ndarray, is_tumor: np.ndarray) -> np.ndarray:
    return (np.log2(np.median(norm[:, is_tumor], axis=1) + 1)
            - np.log2(np.median(norm[:, ~is_tumor], axis=1) + 1))


# ---------------------------------------------------------------- rank tests

def wilcoxon_ranksum(x: np.ndarray, is_tumor: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Two-sided Mann-Whitney U per gene (rows), asymptotic with tie correction.

    No continuity correction, so that a single-stratum van Elteren test gives the
    identical p-value. Returns (pvalue, common-language effect size P(T > N) + 0.5 P(T = N)).
    """
    t, n = x[:, is_tumor], x[:, ~is_tumor]
    res = stats.mannwhitneyu(t, n, axis=1, alternative="two-sided",
                             method="asymptotic", use_continuity=False)
    cles = res.statistic / (t.shape[1] * n.shape[1])
    return res.pvalue, cles


def van_elteren(x: np.ndarray, is_tumor: np.ndarray, strata: np.ndarray) -> dict:
    """van Elteren stratified Wilcoxon test, vectorised over genes (rows of x).

    Within stratum k: W_k = rank sum of tumors, E[W_k] = m_k (n_k + 1) / 2 and Var(W_k)
    is the tie-corrected Wilcoxon variance, written as the exact permutation variance
    m_k (n_k - m_k) / (n_k (n_k - 1)) * sum_i (r_i - (n_k + 1)/2)^2, which equals
    m_k (n_k - m_k) / 12 * [(n_k + 1) - sum_t (t^3 - t) / (n_k (n_k - 1))].
    Weights w_k = 1 / (n_k + 1). Strata with only one condition are skipped.
    """
    x = np.asarray(x, dtype=float)
    is_tumor = np.asarray(is_tumor, dtype=bool)
    strata = np.asarray(strata)
    num = np.zeros(x.shape[0])
    var = np.zeros(x.shape[0])
    n_used = 0
    used = []
    for k in pd.unique(strata):
        idx = strata == k
        m = int(is_tumor[idx].sum())
        n = int(idx.sum())
        if m == 0 or m == n:
            continue
        r = stats.rankdata(x[:, idx], axis=1)
        w_k = r[:, is_tumor[idx]].sum(axis=1)
        e_k = m * (n + 1) / 2
        v_k = m * (n - m) / (n * (n - 1)) * ((r - (n + 1) / 2) ** 2).sum(axis=1)
        wt = 1.0 / (n + 1)
        num += wt * (w_k - e_k)
        var += wt ** 2 * v_k
        n_used += n
        used.append(k)
    with np.errstate(invalid="ignore", divide="ignore"):
        z = num / np.sqrt(var)
    p = 2 * stats.norm.sf(np.abs(z))
    p[var == 0] = np.nan  # gene constant within every informative stratum
    return {"z": z, "pvalue": p, "n_used": n_used, "strata_used": used}


def rank_protocol(norm: pd.DataFrame, is_tumor: np.ndarray, strata=None) -> pd.DataFrame:
    x = norm.to_numpy(dtype=float)
    df = pd.DataFrame(index=norm.index)
    df["log2FC"] = median_log2fc(x, is_tumor)
    t, n = x[:, is_tumor], x[:, ~is_tumor]
    df["CLES"] = stats.mannwhitneyu(t, n, axis=1, method="asymptotic").statistic / (t.shape[1] * n.shape[1])
    if strata is None:
        df["pvalue"], _ = wilcoxon_ranksum(x, is_tumor)
        n_used = x.shape[1]
    else:
        ve = van_elteren(x, is_tumor, strata)
        df["stat_z"] = ve["z"]
        df["pvalue"] = ve["pvalue"]
        n_used = ve["n_used"]
    df = finalize(df)
    df.attrs["n_used"] = n_used
    return df


# ---------------------------------------------------------------- DESeq2

def run_pydeseq2(counts: pd.DataFrame, meta: pd.DataFrame, design: str,
                 n_cpus: int | None = None, shrink: bool = True) -> pd.DataFrame:
    """PyDESeq2 Wald test tumor vs normal; returns MLE and shrunken log2FC.

    Independent filtering is off so FDR is plain BH over all genes, the same as the
    rank protocols. Cook's outlier handling stays at PyDESeq2 defaults.
    """
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats

    dds = DeseqDataSet(counts=counts.T, metadata=meta, design=design,
                       refit_cooks=True, n_cpus=n_cpus, quiet=True)
    dds.deseq2()
    ds = DeseqStats(dds, contrast=["condition", "tumor", "normal"], independent_filter=False,
                    n_cpus=n_cpus, quiet=True)
    ds.summary()
    res = ds.results_df[["baseMean", "log2FoldChange", "lfcSE", "pvalue"]].copy()
    res = res.rename(columns={"log2FoldChange": "log2FC_MLE"})
    if shrink:
        ds.lfc_shrink(coeff=COEF)
        res["log2FC"] = ds.results_df["log2FoldChange"]
    else:
        res["log2FC"] = res["log2FC_MLE"]
    return finalize(res)


def winsorize_counts_global(counts: pd.DataFrame, q: float = 95) -> pd.DataFrame:
    """Old P7 (kept for de_P7_global.csv): cap raw counts at the gene's q-th percentile
    over all samples. Clips the minority condition (normals) on down-regulated genes."""
    cap = np.round(np.percentile(counts.to_numpy(), q, axis=1)).astype(np.int64)
    return pd.DataFrame(np.minimum(counts.to_numpy(), cap[:, None]),
                        index=counts.index, columns=counts.columns)


def winsorize_counts(counts: pd.DataFrame, condition, q: float = 95) -> pd.DataFrame:
    """P7 winsorization following Yang et al. 2024 (Genome Biol 25:282), within condition.

    Counts are divided by DESeq2 size factors; per gene, normalized values above the
    q-th percentile *of the same condition* are replaced by that percentile; values are
    multiplied back by the size factors and rounded. Yang et al. do not state whether the
    percentile is global or per group; per-condition avoids clipping the minority group.
    """
    x = counts.to_numpy(dtype=float)
    sf = size_factors(counts.to_numpy())
    norm = x / sf
    # condition=None: one global percentile over all samples (sensitivity variant "P7g")
    condition = np.zeros(x.shape[1]) if condition is None else np.asarray(condition)
    for c in pd.unique(condition):
        idx = condition == c
        cap = np.percentile(norm[:, idx], q, axis=1)
        norm[:, idx] = np.minimum(norm[:, idx], cap[:, None])
    out = np.round(norm * sf).astype(np.int64)
    return pd.DataFrame(out, index=counts.index, columns=counts.columns)


# ---------------------------------------------------------------- ComBat-seq (R sva)

def find_rscript() -> str:
    for base in [Path(os.environ.get("ProgramFiles", r"C:\Program Files")) / "R",
                 Path(os.environ.get("LOCALAPPDATA", "")) / "Programs" / "R"]:
        hits = sorted(base.glob("R-*/bin/Rscript.exe"))
        if hits:
            return str(hits[-1])
    return "Rscript"


R_LIB = ROOT / ".rlib"  # project-local R package library (sva and dependencies)


def r_env() -> dict:
    """Environment for Rscript: put the project's R library first when it exists, so R
    finds sva no matter which user profile / app sandbox launched the process."""
    env = os.environ.copy()
    if R_LIB.exists():
        env["R_LIBS_USER"] = str(R_LIB)
        env["R_LIBS"] = str(R_LIB)
    return env


COMBAT_R = r'''
suppressMessages(library(sva))
args <- commandArgs(trailingOnly = TRUE)
d <- args[1]
cnt <- as.matrix(read.csv(file.path(d, "counts.csv"), row.names = 1, check.names = FALSE))
meta <- read.csv(file.path(d, "meta.csv"), row.names = 1, check.names = FALSE)
stopifnot(identical(colnames(cnt), rownames(meta)))
sets <- setdiff(colnames(meta), "batch")
for (s in sets) {
  out <- file.path(d, paste0("adjusted_", s, ".csv"))
  if (file.exists(out)) next
  set.seed(42)
  adj <- ComBat_seq(cnt, batch = meta$batch, group = meta[[s]])
  tmp <- paste0(out, ".tmp")
  write.csv(adj, tmp)
  file.rename(tmp, out)
  cat("done", s, "\n")
  flush.console()
}
cat("sva", as.character(packageVersion("sva")), "\n")
'''


def combat_seq_multi(counts: pd.DataFrame, batch: pd.Series, groups: pd.DataFrame,
                     workdir: Path | None = None) -> dict[str, pd.DataFrame]:
    """Run sva::ComBat_seq for several condition label sets in one R session.

    groups: samples x label-sets (column name = set id). Each set is adjusted with
    batch as batch and that set's labels as the protected group covariate. If
    workdir is given, finished sets are kept there and skipped on a rerun, so an
    interrupted run resumes.
    """
    own_tmp = workdir is None
    tmp = Path(tempfile.mkdtemp()) if own_tmp else Path(workdir)
    tmp.mkdir(parents=True, exist_ok=True)
    try:
        if not (tmp / "counts.csv").exists():
            counts.to_csv(tmp / "counts.csv")
        meta = pd.DataFrame({"batch": batch.to_numpy()}, index=counts.columns)
        for c in groups.columns:
            meta[str(c)] = groups[c].to_numpy()
        meta.to_csv(tmp / "meta.csv")
        (tmp / "combat.R").write_text(COMBAT_R, encoding="utf-8")
        out = subprocess.run([find_rscript(), str(tmp / "combat.R"), str(tmp)],
                             capture_output=True, text=True, env=r_env())
        if out.returncode != 0:
            raise RuntimeError(f"ComBat_seq failed:\n{out.stdout}\n{out.stderr}")
        print(out.stdout.strip().splitlines()[-1], flush=True)
        result = {}
        for c in groups.columns:
            adj = pd.read_csv(tmp / f"adjusted_{c}.csv", index_col=0)
            adj.columns = counts.columns
            result[str(c)] = adj.loc[counts.index].round().astype(np.int64)
        return result
    finally:
        if own_tmp:
            import shutil
            shutil.rmtree(tmp, ignore_errors=True)


def combat_seq(counts: pd.DataFrame, batch: pd.Series, condition: pd.Series) -> pd.DataFrame:
    """ComBat-seq for a single label set (batch adjusted, condition kept as group)."""
    groups = pd.DataFrame({"real": condition.to_numpy()}, index=counts.columns)
    return combat_seq_multi(counts, batch, groups)["real"]


# ---------------------------------------------------------------- orchestration

def load_data() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    counts = pd.read_parquet(PROCESSED / "counts_filtered.parquet")
    samples = pd.read_csv(PROCESSED / "samples.csv").set_index("aliquot_barcode")
    samples = samples.loc[counts.columns]
    genes = pd.read_parquet(PROCESSED / "genes.parquet")
    return counts, samples, genes


def meta_for(samples: pd.DataFrame) -> pd.DataFrame:
    meta = samples[["condition", "batch_group"]].copy()
    meta["condition"] = pd.Categorical(meta["condition"], categories=["normal", "tumor"])
    return meta


def get_combat(counts, samples, force=False) -> tuple[pd.DataFrame, float]:
    """ComBat-seq adjusted counts, cached on disk; returns (counts, seconds spent now)."""
    if COMBAT_FILE.exists() and not force:
        timing = json.loads((PROCESSED / "combat_time.json").read_text())
        return pd.read_parquet(COMBAT_FILE), timing["seconds"]
    t0 = time.perf_counter()
    adj = combat_seq(counts, samples["batch_group"], samples["condition"])
    sec = time.perf_counter() - t0
    adj.to_parquet(COMBAT_FILE)
    (PROCESSED / "combat_time.json").write_text(json.dumps({"seconds": sec}))
    return adj, sec


def run_protocol(name: str, counts, samples, n_cpus=None) -> tuple[pd.DataFrame, dict]:
    is_tumor = (samples["condition"] == "tumor").to_numpy()
    meta = meta_for(samples)
    info = {"protocol": name}
    extra = 0.0
    t0 = time.perf_counter()
    if name == "P1":
        res = run_pydeseq2(counts, meta, "~condition", n_cpus)
        info["n_used"] = len(samples)
    elif name == "P2":
        res = run_pydeseq2(counts, meta, "~batch_group + condition", n_cpus)
        info["n_used"] = int(samples["batch_has_both"].sum())
    elif name == "P3":
        adj, extra = get_combat(counts, samples)
        t0 = time.perf_counter()
        res = run_pydeseq2(adj, meta, "~condition", n_cpus)
        info["n_used"] = len(samples)
    elif name == "P4":
        res = rank_protocol(normalize(counts), is_tumor)
        info["n_used"] = res.attrs["n_used"]
    elif name == "P5":
        adj, extra = get_combat(counts, samples)
        t0 = time.perf_counter()
        res = rank_protocol(normalize(adj), is_tumor)
        info["n_used"] = res.attrs["n_used"]
    elif name == "P6":
        res = rank_protocol(normalize(counts), is_tumor, samples["batch_group"].to_numpy())
        info["n_used"] = res.attrs["n_used"]
    elif name == "P7":
        res = run_pydeseq2(winsorize_counts(counts, samples["condition"]), meta,
                           "~condition", n_cpus)
        info["n_used"] = len(samples)
    else:
        raise ValueError(name)
    info["seconds"] = time.perf_counter() - t0 + extra
    info["combat_seconds"] = extra
    info["n_DEG"] = int(res["is_DEG"].sum())
    info["n_up"] = int((res["direction"] == "up").sum())
    info["n_down"] = int((res["direction"] == "down").sum())
    info["n_pvalue_nan"] = int(res["pvalue"].isna().sum())
    return res, info


def save_protocol(name: str, res: pd.DataFrame, genes: pd.DataFrame, info: dict) -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    out = res.copy()
    out.insert(0, "gene_name", genes.loc[out.index, "gene_name"])
    out.index.name = "gene_id"
    out.to_csv(TABLES / f"de_{name}.csv")
    (TABLES / f"de_{name}.json").write_text(json.dumps(info, indent=1))


def load_protocol(name: str) -> pd.DataFrame:
    return pd.read_csv(TABLES / f"de_{name}.csv", index_col="gene_id")


def main(argv=None) -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("protocols", nargs="*", default=PROTOCOLS)
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--n-cpus", type=int, default=max(1, (os.cpu_count() or 2) - 2))
    args = ap.parse_args(argv)
    counts, samples, genes = load_data()
    for name in args.protocols:
        if (TABLES / f"de_{name}.json").exists() and not args.force:
            print(f"{name}: cached, skipping", flush=True)
            continue
        print(f"{name}: running ...", flush=True)
        res, info = run_protocol(name, counts, samples, args.n_cpus)
        save_protocol(name, res, genes, info)
        print(f"{name}: {info}", flush=True)


if __name__ == "__main__":
    main()


# ---------------------------------------------------------------- stage 3 summary

MARKERS = {"COL10A1": "up", "MMP11": "up", "ADH1B": "down", "FABP4": "down", "LEP": "down"}


def summarize_protocols(protocols=PROTOCOLS) -> dict:
    res = {p: load_protocol(p) for p in protocols}
    info = pd.DataFrame([json.loads((TABLES / f"de_{p}.json").read_text()) for p in protocols]) \
        .set_index("protocol")
    deg = {p: set(r.index[r["is_DEG"]]) for p, r in res.items()}
    jac = pd.DataFrame(index=protocols, columns=protocols, dtype=float)
    for a in protocols:
        for b in protocols:
            u = deg[a] | deg[b]
            jac.loc[a, b] = len(deg[a] & deg[b]) / len(u) if u else np.nan
    rows = []
    for p, r in res.items():
        for g, expected in MARKERS.items():
            hit = r[r["gene_name"] == g]
            if hit.empty:
                rows.append({"protocol": p, "gene": g, "expected": expected, "in_filtered_set": False})
                continue
            h = hit.iloc[0]
            rows.append({"protocol": p, "gene": g, "expected": expected, "in_filtered_set": True,
                         "log2FC": h["log2FC"], "FDR": h["FDR"], "is_DEG": bool(h["is_DEG"]),
                         "direction_ok": bool(h["is_DEG"]) and h["direction"] == expected})
    markers = pd.DataFrame(rows)
    info.to_csv(TABLES / "de_protocol_summary.csv")
    jac.to_csv(TABLES / "de_jaccard.csv")
    markers.to_csv(TABLES / "de_markers.csv", index=False)
    return {"info": info, "jaccard": jac, "markers": markers, "results": res}


# ---------------------------------------------------------------- common-LFC sensitivity

# which DESeq2 estimate supplies the common log2FC for each protocol
COMMON_LFC_SOURCE = {"P1": "P1", "P4": "P1", "P7": "P1",   # no batch adjustment
                     "P2": "P2", "P6": "P2",                 # batch in the model / strata
                     "P3": "P3", "P5": "P3"}                 # ComBat-seq adjusted counts
TOP_N = 500


def rank_strength(r: pd.DataFrame) -> pd.Series:
    """Tie-breaker for genes with equal p-value: |Wald z|, |van Elteren z| or |CLES - 0.5|."""
    if "lfcSE" in r:
        return (r["log2FC_MLE"] / r["lfcSE"]).abs()
    if "stat_z" in r:
        return r["stat_z"].abs()
    return (r["CLES"] - 0.5).abs()


def top_genes(r: pd.DataFrame, n: int = TOP_N) -> list:
    order = pd.DataFrame({"p": r["pvalue"].fillna(1.0), "s": -rank_strength(r).fillna(0)})
    return order.sort_values(["p", "s"]).index[:n].tolist()


def jaccard_matrix(sets: dict) -> pd.DataFrame:
    names = list(sets)
    m = pd.DataFrame(index=names, columns=names, dtype=float)
    for a in names:
        for b in names:
            u = sets[a] | sets[b]
            m.loc[a, b] = len(sets[a] & sets[b]) / len(u) if u else np.nan
    return m


def common_lfc_sensitivity(protocols=PROTOCOLS) -> dict:
    out_dir = TABLES / "de_common_lfc"
    out_dir.mkdir(parents=True, exist_ok=True)
    res = {p: load_protocol(p) for p in protocols}
    rows, deg, top = [], {}, {}
    for p in protocols:
        r = res[p]
        lfc = res[COMMON_LFC_SOURCE[p]]["log2FC"].reindex(r.index)
        is_deg = (r["FDR"] < FDR_CUT) & (lfc.abs() > LFC_CUT)
        direction = np.where(is_deg, np.where(lfc > 0, "up", "down"), "")
        pd.DataFrame({"gene_name": r["gene_name"], "common_log2FC": lfc, "FDR": r["FDR"],
                      "is_DEG": is_deg, "direction": direction}).to_csv(out_dir / f"de_{p}_common_lfc.csv")
        deg[p] = set(r.index[is_deg])
        top[p] = set(top_genes(r))
        rows.append({"protocol": p, "lfc_source": COMMON_LFC_SOURCE[p],
                     "n_DEG_own_lfc": int(r["is_DEG"].sum()), "n_DEG_common_lfc": int(is_deg.sum()),
                     "n_up": int((direction == "up").sum()), "n_down": int((direction == "down").sum())})
    summary = pd.DataFrame(rows).set_index("protocol")
    jac = jaccard_matrix(deg)
    jac_top = jaccard_matrix(top)
    overlap_top = pd.DataFrame({a: {b: len(top[a] & top[b]) for b in protocols} for a in protocols})
    summary.to_csv(out_dir / "summary.csv")
    jac.to_csv(out_dir / "jaccard_common_lfc.csv")
    jac_top.to_csv(out_dir / f"jaccard_top{TOP_N}.csv")
    overlap_top.to_csv(out_dir / f"overlap_top{TOP_N}.csv")
    return {"summary": summary, "jaccard": jac, "jaccard_top": jac_top, "overlap_top": overlap_top}
