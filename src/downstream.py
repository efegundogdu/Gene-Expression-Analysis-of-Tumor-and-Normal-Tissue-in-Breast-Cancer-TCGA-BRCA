"""Stage 6: downstream impact of the protocol choice.

1. Enrichment (gseapy, Enrichr libraries GO_Biological_Process_2023 and KEGG_2021_Human,
   run locally against the 15317 tested genes as background) for each protocol's up and
   down genes; Jaccard similarity of the top-20 terms between protocols (Fig 8).
2. Tumor vs normal classification with nested stratified 5-fold CV (L1 logistic
   regression, random forest; class_weight='balanced'). Gene selection is redone inside
   every outer training fold with the protocol's own test, so the test fold never informs
   the gene list (Fig 9).
3. Robust core signature: DEG in P2, P3, P5 and P6 with the same direction, and
   significant in < 5% of B-conf (plate-confounded null) replicates of each of them.
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import de_methods as dm  # noqa: E402

SEED = 42
ROOT = Path(__file__).resolve().parents[1]
PROCESSED = ROOT / "data" / "processed"
TABLES = ROOT / "results" / "tables"
ENRICHR_DIR = PROCESSED / "enrichr"
LIBRARIES = ["GO_Biological_Process_2023", "KEGG_2021_Human"]
TOP_TERMS = 20
ROBUST_PROTOCOLS = ["P2", "P3", "P5", "P6"]
# v2 of the core: B-conf frequency filter without P5, whose own batch-driven false
# positives (ComBat-seq with confounded labels) would otherwise remove genes
ROBUST_FILTER_V2 = ["P2", "P3", "P6"]
ROBUST_MAX_CONF_FREQ = 0.05


# ---------------------------------------------------------------- enrichment

def load_library(name: str) -> dict[str, list[str]]:
    path = ENRICHR_DIR / f"{name}.json"
    if not path.exists():
        import gseapy as gp
        ENRICHR_DIR.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(gp.get_library(name=name, organism="Human")))
    return json.loads(path.read_text())


def enrich(genes: list[str], library: dict, background: list[str]) -> pd.DataFrame:
    """Hypergeometric enrichment (gseapy.enrichr offline mode) sorted by adjusted p."""
    import gseapy as gp
    res = gp.enrichr(gene_list=genes, gene_sets=library, background=background,
                     outdir=None, no_plot=True, verbose=False)
    df = res.results.sort_values(["Adjusted P-value", "P-value"]).reset_index(drop=True)
    # gseapy lists overlap genes in set order, which changes between runs; sort them
    df["Genes"] = df["Genes"].astype(str).str.split(";").apply(lambda g: ";".join(sorted(g)))
    return df


def run_enrichment(protocols=dm.PROTOCOLS) -> dict:
    out_dir = TABLES / "enrichment"
    out_dir.mkdir(parents=True, exist_ok=True)
    genes = pd.read_parquet(PROCESSED / "genes.parquet")
    tested = pd.read_parquet(PROCESSED / "counts_filtered.parquet", columns=[]).index
    background = sorted(set(genes.loc[tested, "gene_name"]))
    libs = {name: load_library(name) for name in LIBRARIES}
    top, n_sig = {}, []
    for p in protocols:
        r = dm.load_protocol(p)
        for direction in ["up", "down"]:
            gl = sorted(set(r.loc[r["direction"] == direction, "gene_name"]))
            for lib_name, lib in libs.items():
                df = enrich(gl, lib, background)
                df.to_csv(out_dir / f"{p}_{direction}_{lib_name}.csv", index=False)
                sig = df[df["Adjusted P-value"] < 0.05]
                top[(lib_name, direction, p)] = list(sig["Term"].head(TOP_TERMS))
                n_sig.append({"protocol": p, "direction": direction, "library": lib_name,
                              "n_genes": len(gl), "n_terms_FDR<0.05": len(sig),
                              "top_term": df["Term"].iloc[0] if len(df) else ""})
    n_sig = pd.DataFrame(n_sig)
    n_sig.to_csv(out_dir / "summary.csv", index=False)
    jac = {}
    for lib_name in LIBRARIES:
        for direction in ["up", "down"]:
            m = pd.DataFrame(index=protocols, columns=protocols, dtype=float)
            for a in protocols:
                for b in protocols:
                    sa, sb = set(top[(lib_name, direction, a)]), set(top[(lib_name, direction, b)])
                    m.loc[a, b] = len(sa & sb) / len(sa | sb) if (sa | sb) else np.nan
            m.to_csv(out_dir / f"jaccard_top{TOP_TERMS}_{lib_name}_{direction}.csv")
            jac[(lib_name, direction)] = m
    with open(out_dir / f"top{TOP_TERMS}_terms.json", "w") as fh:
        json.dump({"|".join(k): v for k, v in top.items()}, fh, indent=1)
    return {"summary": n_sig, "jaccard": jac, "top": top}


# ---------------------------------------------------------------- robust core signature

def robust_core(protocols=ROBUST_PROTOCOLS, max_freq=ROBUST_MAX_CONF_FREQ,
                filter_protocols=None, name: str = "robust_core") -> pd.DataFrame:
    """DEG in every protocol of `protocols` with the same direction, and significant in
    fewer than `max_freq` of the B-conf replicates (Stage 4) for each protocol in
    `filter_protocols` (default: all of `protocols`)."""
    filter_protocols = list(filter_protocols or protocols)
    res = {p: dm.load_protocol(p) for p in protocols}
    bs = pd.read_csv(TABLES / "batch_sensitive_genes.csv", index_col="gene_id")
    genes = res[protocols[0]].index
    dirs = pd.DataFrame({p: res[p]["direction"].fillna("").reindex(genes) for p in protocols})
    same = (dirs.nunique(axis=1) == 1) & (dirs[protocols[0]].isin(["up", "down"]))
    freqs = pd.DataFrame({p: bs.loc[genes, f"freq_conf_{p}"] for p in protocols})
    stable = (freqs[filter_protocols] < max_freq).all(axis=1)
    keep = same & stable
    out = pd.DataFrame({"gene_name": res[protocols[0]].loc[genes, "gene_name"],
                        "direction": dirs[protocols[0]]}, index=genes)
    for p in protocols:
        out[f"log2FC_{p}"] = res[p]["log2FC"]
        out[f"FDR_{p}"] = res[p]["FDR"]
        out[f"freq_conf_{p}"] = freqs[p]
    out = out[keep].sort_values(["direction", "log2FC_P2"], ascending=[False, False])
    out.index.name = "gene_id"
    out.to_csv(TABLES / f"{name}.csv")
    info = {"protocols": protocols, "filter_protocols": filter_protocols, "max_conf_freq": max_freq,
            "DEG_all_same_direction": int(same.sum()),
            "dropped_batch_frequent": int((same & ~stable).sum()),
            "robust_core": int(keep.sum()), "up": int((out["direction"] == "up").sum()),
            "down": int((out["direction"] == "down").sum())}
    (TABLES / f"{name}_info.json").write_text(json.dumps(info, indent=1))
    return out


# ---------------------------------------------------------------- classification

CLASSIF_DIR = PROCESSED / "classif"
SELECTORS = ["P1", "P2", "P4", "P6", "P7"]  # P3/P5 need ComBat-seq per fold (see README)
N_OUTER, N_INNER = 5, 5


def train_size_factors(train_counts: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Median-of-ratios reference (log geometric means) from training samples only."""
    with np.errstate(divide="ignore"):
        lg = np.log(train_counts.astype(float)).mean(axis=1)
    return lg, np.isfinite(lg)


def apply_size_factors(counts: np.ndarray, lg: np.ndarray, ok: np.ndarray) -> np.ndarray:
    with np.errstate(divide="ignore"):
        lc = np.log(counts.astype(float))
    ratios = lc[ok] - lg[ok, None]
    ratios[~np.isfinite(ratios)] = np.nan  # zero counts in a sample are ignored
    return np.exp(np.nanmedian(ratios, axis=0))


def select_genes(protocol: str, counts: pd.DataFrame, samples: pd.DataFrame,
                 n_cpus: int | None = None) -> pd.Series:
    """DEG direction ('up'/'down') on the given (training) samples with the protocol's test."""
    is_tumor = (samples["condition"] == "tumor").to_numpy()
    meta = dm.meta_for(samples)
    if protocol == "P1":
        res = dm.run_pydeseq2(counts, meta, "~condition", n_cpus)
    elif protocol == "P2":
        res = dm.run_pydeseq2(counts, meta, "~batch_group + condition", n_cpus)
    elif protocol == "P4":
        res = dm.rank_protocol(dm.normalize(counts), is_tumor)
    elif protocol == "P6":
        res = dm.rank_protocol(dm.normalize(counts), is_tumor, samples["batch_group"].to_numpy())
    elif protocol == "P7":
        res = dm.run_pydeseq2(dm.winsorize_counts(counts, samples["condition"]), meta,
                              "~condition", n_cpus)
    else:
        raise ValueError(protocol)
    return res.loc[res["is_DEG"], "direction"]


def make_models():
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.linear_model import LogisticRegression
    from sklearn.model_selection import GridSearchCV, StratifiedKFold
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler
    inner = StratifiedKFold(N_INNER, shuffle=True, random_state=SEED)
    lr = GridSearchCV(make_pipeline(StandardScaler(), LogisticRegression(
        penalty="l1", solver="liblinear", class_weight="balanced", max_iter=5000,
        random_state=SEED)),
        {"logisticregression__C": [0.01, 0.1, 1.0]}, cv=inner, scoring="roc_auc", n_jobs=-1)
    rf = GridSearchCV(RandomForestClassifier(n_estimators=300, class_weight="balanced",
                                             random_state=SEED, n_jobs=-1),
                      {"max_features": ["sqrt", 0.02]}, cv=inner, scoring="roc_auc", n_jobs=1)
    return {"L1-LR": lr, "RF": rf}


def evaluate_fold(model, X_tr, y_tr, X_te, y_te, feature_names) -> dict:
    from sklearn.metrics import f1_score, roc_auc_score
    model.fit(X_tr, y_tr)
    prob = model.predict_proba(X_te)[:, 1]
    pred = (prob >= 0.5).astype(int)
    est = model.best_estimator_
    if hasattr(est, "named_steps"):
        coef = est.named_steps["logisticregression"].coef_.ravel()
        used = [g for g, c in zip(feature_names, coef) if c != 0]
    else:
        imp = est.feature_importances_
        used = [feature_names[i] for i in np.argsort(imp)[::-1][:50]]
    params = {k: (v if isinstance(v, (int, float, str)) else str(v))
              for k, v in model.best_params_.items()}
    return {"roc_auc": float(roc_auc_score(y_te, prob)),
            "f1_tumor": float(f1_score(y_te, pred, pos_label=1)),
            "f1_normal": float(f1_score(y_te, pred, pos_label=0)),
            "best_params": params, "model_genes": used}


def log_features(counts: np.ndarray, sf: np.ndarray) -> np.ndarray:
    return np.log2(counts / sf + 1).T  # samples x genes


def run_classification(selectors=SELECTORS, fixed_lists: dict | None = None,
                       n_cpus: int = 12, log=print) -> pd.DataFrame:
    """Nested CV. selectors: gene selection inside each outer training fold.
    fixed_lists: name -> gene list chosen beforehand on all data (not nested; see README)."""
    from sklearn.model_selection import StratifiedKFold
    counts, samples, _ = dm.load_data()
    y = (samples["condition"] == "tumor").astype(int).to_numpy()
    outer = StratifiedKFold(N_OUTER, shuffle=True, random_state=SEED)
    CLASSIF_DIR.mkdir(parents=True, exist_ok=True)
    jobs = [(p, "nested") for p in selectors] + [(n, "fixed") for n in (fixed_lists or {})]
    for k, (tr, te) in enumerate(outer.split(np.zeros(len(y)), y)):
        c_tr, c_te = counts.iloc[:, tr], counts.iloc[:, te]
        lg, ok = train_size_factors(c_tr.to_numpy())
        sf_tr = apply_size_factors(c_tr.to_numpy(), lg, ok)
        sf_te = apply_size_factors(c_te.to_numpy(), lg, ok)
        for name, kind in jobs:
            path = CLASSIF_DIR / f"{name}_fold{k}.json"
            if path.exists():
                continue
            t0 = time.perf_counter()
            if kind == "nested":
                genes = list(select_genes(name, c_tr, samples.iloc[tr], n_cpus).index)
            else:
                genes = [g for g in fixed_lists[name] if g in counts.index]
            t_sel = time.perf_counter() - t0
            rows = counts.index.get_indexer(genes)
            X_tr = log_features(c_tr.to_numpy()[rows], sf_tr)
            X_te = log_features(c_te.to_numpy()[rows], sf_te)
            rec = {"name": name, "kind": kind, "fold": k, "n_genes": len(genes),
                   "selected_genes": genes, "seconds_selection": t_sel, "models": {}}
            for mname, model in make_models().items():
                rec["models"][mname] = evaluate_fold(model, X_tr, y[tr], X_te, y[te], genes)
            rec["seconds_total"] = time.perf_counter() - t0
            path.write_text(json.dumps(rec))
            log(f"fold {k} {name}: {len(genes)} genes, "
                f"AUC LR {rec['models']['L1-LR']['roc_auc']:.4f} "
                f"RF {rec['models']['RF']['roc_auc']:.4f}, {rec['seconds_total']:.0f}s")
    return summarize_classification()


def _pairwise_jaccard(sets: list[set]) -> float:
    vals = [len(a & b) / len(a | b) for i, a in enumerate(sets) for b in sets[i + 1:] if a | b]
    return float(np.mean(vals)) if vals else np.nan


def summarize_classification() -> pd.DataFrame:
    recs = [json.loads(p.read_text()) for p in sorted(CLASSIF_DIR.glob("*_fold*.json"))]
    rows, stab = [], []
    for r in recs:
        for m, v in r["models"].items():
            rows.append({"name": r["name"], "kind": r["kind"], "fold": r["fold"], "model": m,
                         "n_genes": r["n_genes"], "n_model_genes": len(v["model_genes"]),
                         "roc_auc": v["roc_auc"], "f1_tumor": v["f1_tumor"],
                         "f1_normal": v["f1_normal"], "best_params": json.dumps(v["best_params"])})
    per_fold = pd.DataFrame(rows)
    per_fold.to_csv(TABLES / "classification_per_fold.csv", index=False)
    for name in sorted({r["name"] for r in recs}):
        rr = [r for r in recs if r["name"] == name]
        row = {"name": name,
               "stability_selected": _pairwise_jaccard([set(r["selected_genes"]) for r in rr])}
        for m in ["L1-LR", "RF"]:
            row[f"stability_{m}"] = _pairwise_jaccard([set(r["models"][m]["model_genes"]) for r in rr])
        stab.append(row)
    stab = pd.DataFrame(stab).set_index("name")
    summ = per_fold.groupby(["name", "model"]).agg(
        kind=("kind", "first"), folds=("fold", "count"), n_genes=("n_genes", "mean"),
        n_model_genes=("n_model_genes", "mean"), roc_auc=("roc_auc", "mean"),
        roc_auc_sd=("roc_auc", "std"), f1_tumor=("f1_tumor", "mean"),
        f1_normal=("f1_normal", "mean"), f1_normal_sd=("f1_normal", "std")).reset_index()
    summ = summ.join(stab, on="name")
    summ.to_csv(TABLES / "classification_summary.csv", index=False)
    return summ


def run_all(n_cpus: int = 12) -> dict:
    """Stage 6 end to end; classification folds already on disk are reused."""
    import plots
    enr = run_enrichment()
    fig8 = plots.fig8_term_jaccard(enr["jaccard"], TOP_TERMS)
    core = robust_core()
    core_v2 = robust_core(filter_protocols=ROBUST_FILTER_V2, name="robust_core_v2")
    p1 = dm.load_protocol("P1")
    fixed = {"core_fixed": core.index.tolist(), "core_v2_fixed": core_v2.index.tolist(),
             "P1_fixed": p1.index[p1["is_DEG"]].tolist()}
    summ = run_classification(fixed_lists=fixed, n_cpus=n_cpus)
    per_fold = pd.read_csv(TABLES / "classification_per_fold.csv")
    fig9 = plots.fig9_classification(per_fold, summ)
    return {"enrichment": enr, "fig8": fig8, "robust_core": core, "robust_core_v2": core_v2,
            "classification": summ,
            "fig9": fig9}


if __name__ == "__main__":
    import warnings
    warnings.simplefilter("ignore")
    out = run_all()
    print(json.loads((TABLES / "robust_core_info.json").read_text()))
    print(out["classification"].round(4).to_string(index=False))
