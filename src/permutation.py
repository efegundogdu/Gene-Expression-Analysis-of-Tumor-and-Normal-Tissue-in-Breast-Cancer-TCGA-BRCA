"""Stage 4: null datasets and false-positive counting.

Null A  (method-driven, Li et al. replication): real samples, condition labels shuffled
        over all samples.
Null B  (batch-driven, main novel analysis): tumors only (one per patient), so every DEG
        is a false positive. Fake "normal" labels are given
          B-rand: at random, independent of plate;
          B-conf: mimicking the real plate x condition structure (~20% of fake normals on
                  plates that become fake-normal-only, the rest spread over "mixed" plates
                  with real normal fractions).
        FP_batch = FP_conf - FP_rand.
Null C  (secondary): restricted permutation within batch groups, using only samples in the
        15 batch groups that hold both conditions.

Rank protocols (P4, P6) run on the full sets; DESeq2 protocols (P1, P2, P3, P7) and P5
(needs ComBat-seq) run on balanced subsamples. P3 and P5 share one ComBat-seq run per
replicate. Per-gene p-values and log2FC are stored per replicate so a run can resume.
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
import batch as batchmod  # noqa: E402
import de_methods as dm  # noqa: E402

SEED = 42
ROOT = Path(__file__).resolve().parents[1]
PROCESSED = ROOT / "data" / "processed"
PERM = PROCESSED / "perm"
TABLES = ROOT / "results" / "tables"
LOGS = ROOT / "results" / "logs"

SUB_TUMOR_RATIO = 3          # null A subsample: all normals + 3x tumors
SUB_B_SIZE = 452             # null B subsample: tumors only
FAKE_NORMAL_FRAC = 113 / 1208  # same normal fraction as the real data
V_TOL = 0.03                 # B-conf: |V - V_real| tolerance
DESEQ_SET = ["P1", "P2", "P3", "P5", "P7"]
RANK_SET = ["P4", "P6"]


# ---------------------------------------------------------------- sample sets

def _streams(n: int):
    """Independent, reproducible RNG streams derived from seed 42."""
    return [np.random.default_rng(s) for s in np.random.SeedSequence(SEED).spawn(n)]


def build_sample_sets(samples: pd.DataFrame) -> dict[str, pd.DataFrame]:
    """Return named sample tables (index = aliquot barcode)."""
    rng_a, rng_b = _streams(2)
    s = samples.sort_index()
    tumors = s[s["condition"] == "tumor"]
    normals = s[s["condition"] == "normal"]
    sub_a = pd.concat([normals, tumors.loc[rng_a.choice(tumors.index, SUB_TUMOR_RATIO * len(normals),
                                                        replace=False)]]).sort_index()
    sub_b = tumors.loc[rng_b.choice(tumors.index, SUB_B_SIZE, replace=False)].sort_index()
    mixed_groups = set(s.loc[s["batch_has_both"], "batch_group"])
    sets = {
        "full": s,
        "subA": sub_a,
        "tumors": tumors.copy(),
        "subB": sub_b.copy(),
        "mixed_full": s[s["batch_group"].isin(mixed_groups)],
        "mixed_subA": sub_a[sub_a["batch_group"].isin(mixed_groups)],
    }
    # null B: batch groups rebuilt with the same rule on the tumor-only sets
    for name in ["tumors", "subB"]:
        mapping = batchmod.group_plates(sets[name])
        sets[name]["batch_group"] = sets[name]["plate"].map(mapping)
    return sets


# ---------------------------------------------------------------- label generators

def cramers_v(groups, labels) -> float:
    tab = pd.crosstab(np.asarray(groups), np.asarray(labels))
    if tab.shape[1] < 2:
        return 0.0
    chi2 = stats.chi2_contingency(tab, correction=False)[0]
    return float(np.sqrt(chi2 / (tab.to_numpy().sum() * (min(tab.shape) - 1))))


def real_plate_structure(samples: pd.DataFrame) -> dict:
    tab = pd.crosstab(samples["plate"], samples["condition"])
    n_normal = int(tab["normal"].sum())
    pure = tab[(tab["tumor"] == 0)]
    mixed = tab[(tab["tumor"] > 0) & (tab["normal"] > 0)]
    tumor_plates = tab[tab["tumor"] > 0]
    return {"frac_pure": float(pure["normal"].sum() / n_normal),
            "mixed_fracs": (mixed["normal"] / mixed.sum(axis=1)).to_numpy(),
            "frac_plates_mixed": len(mixed) / len(tumor_plates),
            "V_real": cramers_v(samples["plate"], samples["condition"])}


def null_a(df: pd.DataFrame, rng) -> np.ndarray:
    return rng.permutation((df["condition"] == "tumor").to_numpy())


def null_c(df: pd.DataFrame, rng) -> np.ndarray:
    lab = (df["condition"] == "tumor").to_numpy().copy()
    for _, idx in df.groupby("batch_group").indices.items():
        lab[idx] = rng.permutation(lab[idx])
    return lab


def null_b_rand(df: pd.DataFrame, n_fake: int, rng) -> np.ndarray:
    is_tumor = np.ones(len(df), bool)
    is_tumor[rng.choice(len(df), n_fake, replace=False)] = False
    return is_tumor


def null_b_conf(df: pd.DataFrame, n_fake: int, ref: dict, rng) -> tuple[np.ndarray, dict]:
    """Fake normals placed to mimic the real plate x condition structure."""
    plates = df["plate"].to_numpy()
    sizes = pd.Series(plates).value_counts()
    is_tumor = np.ones(len(df), bool)
    # 1) whole small plates become fake-normal-only (~frac_pure of fake normals)
    target_pure = int(round(ref["frac_pure"] * n_fake))
    pure, remaining = [], target_pure
    for p in rng.permutation(sizes.index.to_numpy()):
        if remaining == 0:
            break
        if sizes[p] <= remaining:
            pure.append(p)
            remaining -= sizes[p]
    for p in pure:
        is_tumor[plates == p] = False
    n_pure = int((~is_tumor).sum())
    # 2) the rest goes to "mixed" plates with real normal fractions
    rest = n_fake - n_pure
    others = np.array([p for p in sizes.index if p not in pure and sizes[p] >= 2])
    k = max(1, int(round(ref["frac_plates_mixed"] * len(others))))
    chosen = rng.choice(others, size=min(k, len(others)), replace=False)
    fr = rng.choice(ref["mixed_fracs"], size=len(chosen), replace=len(chosen) > len(ref["mixed_fracs"]))
    cap = np.array([sizes[p] - 1 for p in chosen])  # keep >= 1 tumor so the plate stays mixed
    want = fr * np.array([sizes[p] for p in chosen])
    want = want / want.sum() * rest
    alloc = np.minimum(np.floor(want).astype(int), cap)
    while alloc.sum() < rest:
        room = alloc < cap
        if not room.any():
            break
        gap = np.where(room, want - alloc, -np.inf)
        alloc[int(np.argmax(gap))] += 1
    for p, a in zip(chosen, alloc):
        idx = np.flatnonzero(plates == p)
        is_tumor[rng.choice(idx, a, replace=False)] = False
    info = {"n_fake": int((~is_tumor).sum()), "n_pure_plate": n_pure, "pure_plates": list(pure),
            "n_mixed_plates": int((alloc > 0).sum())}
    return is_tumor, info


def generate_labels(sets: dict, samples_all: pd.DataFrame, B: int) -> pd.DataFrame:
    """Draw and save all label sets; returns per-replicate diagnostics."""
    PERM.mkdir(parents=True, exist_ok=True)
    ref = real_plate_structure(samples_all)
    rng_a_full, rng_a_sub, rng_b_full, rng_b_sub, rng_c_full, rng_c_sub = _streams(8)[2:]
    plans = {
        "A_full": ("full", lambda df, r: null_a(df, r), rng_a_full),
        "A_sub": ("subA", lambda df, r: null_a(df, r), rng_a_sub),
        "C_full": ("mixed_full", lambda df, r: null_c(df, r), rng_c_full),
        "C_sub": ("mixed_subA", lambda df, r: null_c(df, r), rng_c_sub),
    }
    diag = []
    for name, (set_name, fn, rng) in plans.items():
        df = sets[set_name]
        truth = (df["condition"] == "tumor").to_numpy()
        mat = np.column_stack([fn(df, rng) for _ in range(B)])
        pd.DataFrame(mat, index=df.index, columns=[f"r{b:03d}" for b in range(B)]).to_parquet(
            PERM / f"labels_{name}.parquet")
        for b in range(B):
            diag.append({"null": name, "rep": b, "n": len(df), "n_fake_normal": int((~mat[:, b]).sum()),
                         "corr_with_truth": float(np.corrcoef(mat[:, b], truth)[0, 1]),
                         "V_plate": cramers_v(df["plate"], mat[:, b])})
    # null B: the same sample set and the same number of fake normals in rand and conf
    # target V for B-conf: real-label V on the real set of matching size
    v_target = {"tumors": cramers_v(sets["full"]["plate"], sets["full"]["condition"]),
                "subB": cramers_v(sets["subA"]["plate"], sets["subA"]["condition"])}
    for set_name, rng, tag in [("tumors", rng_b_full, "full"), ("subB", rng_b_sub, "sub")]:
        df = sets[set_name]
        n_fake = int(round(FAKE_NORMAL_FRAC * len(df)))
        conf_cols, rand_cols = [], []
        for b in range(B):
            # rejection sampling: keep a B-conf draw only if its plate x label V is
            # within V_TOL of the real data's V
            for tries in range(1, 5001):
                conf, info = null_b_conf(df, n_fake, ref, rng)
                if abs(cramers_v(df["plate"], conf) - v_target[set_name]) <= V_TOL:
                    break
            else:
                raise RuntimeError(f"B-conf: no draw within V tolerance for {set_name}")
            # rand uses exactly as many fake normals as conf produced
            rand = null_b_rand(df, info["n_fake"], rng)
            conf_cols.append(conf)
            rand_cols.append(rand)
            for nm, lab in [(f"Bconf_{tag}", conf), (f"Brand_{tag}", rand)]:
                groups = df["batch_group"].to_numpy()
                single = pd.crosstab(groups, lab)
                diag.append({"null": nm, "rep": b, "n": len(df), "n_fake_normal": int((~lab).sum()),
                             "V_plate": cramers_v(df["plate"], lab),
                             "n_pure_plate_fake": info["n_pure_plate"] if nm.startswith("Bconf") else np.nan,
                             "V_target": v_target[set_name], "draws": tries if nm.startswith("Bconf") else np.nan,
                             "n_batch_groups": int(single.shape[0]),
                             "n_single_condition_groups": int((single.min(axis=1) == 0).sum())
                             if single.shape[1] == 2 else int(single.shape[0])})
        for nm, cols in [(f"Bconf_{tag}", conf_cols), (f"Brand_{tag}", rand_cols)]:
            pd.DataFrame(np.column_stack(cols), index=df.index,
                         columns=[f"r{b:03d}" for b in range(B)]).to_parquet(PERM / f"labels_{nm}.parquet")
    diag = pd.DataFrame(diag)
    diag.to_csv(PERM / "label_diagnostics.csv", index=False)
    (PERM / "reference.json").write_text(json.dumps(
        {k: (v.tolist() if isinstance(v, np.ndarray) else v) for k, v in ref.items()}, indent=1))
    return diag


# ---------------------------------------------------------------- datasets

# dataset name -> (sample set, protocols)
DATASETS = {
    "A_full": ("full", RANK_SET),
    "A_sub": ("subA", DESEQ_SET + RANK_SET),
    "Brand_full": ("tumors", RANK_SET),
    "Bconf_full": ("tumors", RANK_SET),
    "Brand_sub": ("subB", DESEQ_SET + RANK_SET),
    "Bconf_sub": ("subB", DESEQ_SET + RANK_SET),
    "C_full": ("mixed_full", RANK_SET),
    "C_sub": ("mixed_subA", DESEQ_SET + RANK_SET),
}
# real-label analyses on the same sample sets (for FDR estimates and the subsample check)
REAL_DATASETS = {"real_subA": ("subA", DESEQ_SET + RANK_SET),
                 "real_mixed_full": ("mixed_full", RANK_SET),
                 "real_mixed_subA": ("mixed_subA", DESEQ_SET + RANK_SET)}


def result_path(dataset: str, protocol: str, rep: str) -> Path:
    return PERM / "results" / dataset / protocol / f"{rep}.npz"


def save_result(path: Path, res: pd.DataFrame, seconds: float) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(".tmp.npz")
    np.savez(tmp, pvalue=res["pvalue"].to_numpy(np.float64),
             log2FC=res["log2FC"].to_numpy(np.float32), seconds=np.float64(seconds))
    os.replace(tmp, path)


def run_one(protocol: str, counts: pd.DataFrame, df: pd.DataFrame, is_tumor: np.ndarray,
            norm: pd.DataFrame | None = None, adj: pd.DataFrame | None = None,
            n_cpus: int | None = None) -> pd.DataFrame:
    """One protocol on one label set. Returns log2FC / pvalue / FDR / is_DEG."""
    meta = pd.DataFrame({"condition": pd.Categorical(np.where(is_tumor, "tumor", "normal"),
                                                     categories=["normal", "tumor"]),
                         "batch_group": df["batch_group"].to_numpy()}, index=counts.columns)
    if protocol == "P1":
        return dm.run_pydeseq2(counts, meta, "~condition", n_cpus)
    if protocol == "P2":
        return dm.run_pydeseq2(counts, meta, "~batch_group + condition", n_cpus)
    if protocol == "P3":
        return dm.run_pydeseq2(adj, meta, "~condition", n_cpus)
    if protocol == "P4":
        return dm.rank_protocol(norm, is_tumor)
    if protocol == "P5":
        return dm.rank_protocol(dm.normalize(adj), is_tumor)
    if protocol == "P6":
        return dm.rank_protocol(norm, is_tumor, df["batch_group"].to_numpy())
    if protocol == "P7":
        return dm.run_pydeseq2(dm.winsorize_counts(counts, meta["condition"].to_numpy()), meta,
                               "~condition", n_cpus)
    raise ValueError(protocol)


def load_context():
    counts_all = pd.read_parquet(PROCESSED / "counts_filtered.parquet")
    samples = pd.read_csv(PROCESSED / "samples.csv").set_index("aliquot_barcode")
    sets = build_sample_sets(samples)
    return counts_all, samples, sets


def labels_for(dataset: str, sets: dict) -> pd.DataFrame:
    if dataset.startswith("real_"):
        df = sets[REAL_DATASETS[dataset][0]]
        return pd.DataFrame({"real": (df["condition"] == "tumor").to_numpy()}, index=df.index)
    return pd.read_parquet(PERM / f"labels_{dataset}.parquet")


def all_datasets():
    return {**{k: v for k, v in REAL_DATASETS.items()}, **DATASETS}


# ---------------------------------------------------------------- ComBat-seq workers

def free_ram_gb() -> tuple[float, float]:
    """(total, available) physical memory in GB (Windows via ctypes, else /proc)."""
    if os.name == "nt":
        import ctypes

        class MS(ctypes.Structure):
            _fields_ = [("dwLength", ctypes.c_ulong), ("dwMemoryLoad", ctypes.c_ulong),
                        ("ullTotalPhys", ctypes.c_ulonglong), ("ullAvailPhys", ctypes.c_ulonglong),
                        ("ullTotalPageFile", ctypes.c_ulonglong), ("ullAvailPageFile", ctypes.c_ulonglong),
                        ("ullTotalVirtual", ctypes.c_ulonglong), ("ullAvailVirtual", ctypes.c_ulonglong),
                        ("ullAvailExtendedVirtual", ctypes.c_ulonglong)]
        m = MS()
        m.dwLength = ctypes.sizeof(MS)
        ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(m))
        return m.ullTotalPhys / 2**30, m.ullAvailPhys / 2**30
    info = dict(line.split(":") for line in open("/proc/meminfo"))
    return int(info["MemTotal"].split()[0]) / 2**20, int(info["MemAvailable"].split()[0]) / 2**20


def n_workers_for_ram(per_worker_gb: float, reserve_gb: float, cap: int, frac: float = 0.70) -> int:
    """Workers so that total RAM use stays <= frac of physical RAM."""
    total, avail = free_ram_gb()
    budget = frac * total - (total - avail) - reserve_gb
    return int(max(1, min(cap, budget // per_worker_gb)))


def kill_tree(proc: subprocess.Popen) -> None:
    """Stop an R worker and its children (Rscript.exe starts a child R process on Windows)."""
    if os.name == "nt":
        subprocess.run(["taskkill", "/T", "/F", "/PID", str(proc.pid)], capture_output=True)
    else:
        proc.terminate()
    try:
        proc.wait(timeout=30)
    except subprocess.TimeoutExpired:
        proc.kill()


def launch_combat_worker(workdir: Path, counts: pd.DataFrame, batch: np.ndarray,
                         label_sets: pd.DataFrame) -> subprocess.Popen:
    workdir.mkdir(parents=True, exist_ok=True)
    if not (workdir / "counts.csv").exists():
        counts.to_csv(workdir / "counts.csv")
    meta = pd.DataFrame({"batch": batch}, index=counts.columns)
    for c in label_sets.columns:
        meta[c] = np.where(label_sets[c].to_numpy(), "tumor", "normal")
    meta.to_csv(workdir / "meta.csv")
    (workdir / "combat.R").write_text(dm.COMBAT_R, encoding="utf-8")
    log = open(workdir / "worker.log", "a")
    return subprocess.Popen([dm.find_rscript(), str(workdir / "combat.R"), str(workdir)],
                            stdout=log, stderr=subprocess.STDOUT, env=dm.r_env())


# ---------------------------------------------------------------- main runner

STOP_FILE = PERM / "STOP"
RUN_CONFIG = PERM / "run_config.json"


def target_reps(ds: str, b_full: int, b_sub: int) -> list[str]:
    """Replicates to run: rank-only full-data sets use b_full, subsample sets use b_sub."""
    if ds.startswith("real_"):
        return ["real"]
    b = b_full if ds.endswith("_full") else b_sub
    return [f"r{i:03d}" for i in range(b)]


def pending_jobs(datasets: list[str], sets: dict, b_full: int, b_sub: int) -> list[tuple]:
    """(dataset, protocol, rep) triples without a saved result."""
    jobs = []
    for ds in datasets:
        set_name, protocols = all_datasets()[ds]
        cols = [c for c in target_reps(ds, b_full, b_sub) if c in labels_for(ds, sets).columns]
        for rep in cols:
            for p in protocols:
                if not result_path(ds, p, rep).exists():
                    jobs.append((ds, p, rep))
    return jobs


def run(datasets: list[str], b_full: int, b_sub: int, n_cpus: int, max_workers: int,
        r_worker_gb: float, reserve_gb: float, log=print) -> None:
    counts_all, samples, sets = load_context()
    if STOP_FILE.exists():
        STOP_FILE.unlink()
    # unprocessed ComBat-seq outputs from an interrupted run are simply recomputed
    shutil.rmtree(PERM / "combat_work", ignore_errors=True)
    cache = {}

    def ctx(ds):
        if ds not in cache:
            set_name = all_datasets()[ds][0]
            df = sets[set_name]
            counts = counts_all[df.index]
            cache[ds] = (df, counts, dm.normalize(counts), labels_for(ds, sets))
        return cache[ds]

    jobs = pending_jobs(datasets, sets, b_full, b_sub)
    combat_jobs = sorted({(d, r) for d, p, r in jobs if p in ("P3", "P5")})
    other_jobs = [j for j in jobs if j[1] not in ("P3", "P5")]
    log(f"pending: {len(jobs)} protocol runs; ComBat-seq label sets: {len(combat_jobs)}")

    # start R workers: one chunk of label sets per worker, grouped by dataset
    workers = []
    if combat_jobs:
        k = n_workers_for_ram(r_worker_gb, reserve_gb=reserve_gb, cap=max_workers)
        log(f"ComBat-seq workers: {k} (RAM total/avail = {free_ram_gb()[0]:.1f}/{free_ram_gb()[1]:.1f} GB)")
        chunks = [combat_jobs[i::k] for i in range(k)]
        for wi, chunk in enumerate(chunks):
            by_ds = {}
            for d, r in chunk:
                by_ds.setdefault(d, []).append(r)
            workers.append({"queue": list(by_ds.items()), "proc": None, "dir": None, "id": wi})
    combat_dir = PERM / "combat_work"

    def start_next(w):
        while w["queue"]:
            d, reps_d = w["queue"].pop(0)
            df, counts, _, labels = ctx(d)
            wd = combat_dir / f"w{w['id']}_{d}"
            w["proc"] = launch_combat_worker(wd, counts, df["batch_group"].to_numpy(), labels[reps_d])
            w["dir"], w["ds"] = wd, d
            return
        w["proc"] = None

    for w in workers:
        start_next(w)

    def process_combat_outputs():
        did = False
        for w in workers:
            if w["dir"] is None:
                continue
            for f in sorted(w["dir"].glob("adjusted_*.csv")):
                rep = f.stem.replace("adjusted_", "")
                d = w["ds"]
                df, counts, _, labels = ctx(d)
                adj = pd.read_csv(f, index_col=0)
                adj.columns = counts.columns
                adj = adj.loc[counts.index].round().astype(np.int64)
                for p in ("P5", "P3"):
                    if p in all_datasets()[d][1] and not result_path(d, p, rep).exists():
                        t0 = time.perf_counter()
                        res = run_one(p, counts, df, labels[rep].to_numpy(), adj=adj, n_cpus=n_cpus)
                        save_result(result_path(d, p, rep), res, time.perf_counter() - t0)
                        log(f"{d} {p} {rep} FP={int((res['FDR'] < 0.05).sum())} "
                            f"t={time.perf_counter() - t0:.0f}s")
                f.unlink()
                did = True
        return did

    while other_jobs or any(w["proc"] is not None for w in workers):
        if STOP_FILE.exists():
            for w in workers:
                if w["proc"] is not None and w["proc"].poll() is None:
                    kill_tree(w["proc"])
            log("STOP file found: R workers terminated, exiting. Rerun the same command to resume.")
            return
        if process_combat_outputs():
            continue
        for w in workers:
            if w["proc"] is not None and w["proc"].poll() is not None:
                rc = w["proc"].returncode
                process_combat_outputs()
                if rc != 0:
                    log(f"R worker {w['id']} exited with {rc}; see {w['dir'] / 'worker.log'}")
                    for other in workers:
                        if other["proc"] is not None and other["proc"].poll() is None:
                            kill_tree(other["proc"])
                    log("ERROR: ComBat-seq failed, stopping the run (fix the R error, then resume).")
                    raise SystemExit(2)
                start_next(w)
        if other_jobs:
            d, p, rep = other_jobs.pop(0)
            df, counts, norm, labels = ctx(d)
            t0 = time.perf_counter()
            res = run_one(p, counts, df, labels[rep].to_numpy(), norm=norm, n_cpus=n_cpus)
            save_result(result_path(d, p, rep), res, time.perf_counter() - t0)
            log(f"{d} {p} {rep} FP={int((res['FDR'] < 0.05).sum())} t={time.perf_counter() - t0:.1f}s")
        else:
            time.sleep(2)
    process_combat_outputs()
    log("all done")


def print_status() -> None:
    """Completed replicates per dataset and protocol vs the target of the last run."""
    cfg = json.loads(RUN_CONFIG.read_text()) if RUN_CONFIG.exists() else {"b_full": 100, "b_sub": 50}
    rows = []
    for ds, (_, protocols) in all_datasets().items():
        target = len(target_reps(ds, cfg["b_full"], cfg["b_sub"]))
        for p in protocols:
            d = PERM / "results" / ds / p
            done = len(list(d.glob("*.npz"))) if d.exists() else 0
            rows.append({"dataset": ds, "protocol": p, "done": done, "target": target})
    df = pd.DataFrame(rows)
    tab = df.pivot(index="dataset", columns="protocol", values="done").astype("Int64")
    print(tab.astype(str).replace("<NA>", "-").to_string())
    print()
    print(f"total: {int(df['done'].sum())} / {int(df['target'].sum())} protocol runs "
          f"(b_full={cfg['b_full']}, b_sub={cfg['b_sub']})")
    print("complete" if (df["done"] >= df["target"]).all() else "incomplete")


def main(argv=None):
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    g = sub.add_parser("labels")
    g.add_argument("--B", type=int, default=100)
    r = sub.add_parser("run")
    r.add_argument("datasets", nargs="*", default=list(REAL_DATASETS) + list(DATASETS))
    r.add_argument("--b-full", type=int, default=100, help="replicates for full-data rank sets")
    r.add_argument("--b-sub", type=int, default=50, help="replicates for subsample sets")
    r.add_argument("--n-cpus", type=int, default=8)
    r.add_argument("--max-workers", type=int, default=6)
    r.add_argument("--r-worker-gb", type=float, default=1.3)
    r.add_argument("--reserve-gb", type=float, default=2.0)
    st = sub.add_parser("status")
    sm = sub.add_parser("summary")
    args = ap.parse_args(argv)
    if args.cmd == "status":
        print_status()
        return
    if args.cmd == "summary":
        res = run_summary()
        print(res["summary"][["protocol", "threshold", "FP_A", "FP_rand", "FP_conf", "FP_batch",
                              "FP_C", "FDR_A"]].round(3).to_string(index=False))
        return
    if args.cmd == "labels":
        _, samples, sets = load_context()
        diag = generate_labels(sets, samples, args.B)
        print(diag.groupby("null").agg(n=("n", "first"), fake_normal=("n_fake_normal", "mean"),
                                       V=("V_plate", "mean"), V_sd=("V_plate", "std"),
                                       corr=("corr_with_truth", "mean")).round(3).to_string())
    else:
        LOGS.mkdir(parents=True, exist_ok=True)
        logf = open(LOGS / "permutation.log", "a", encoding="utf-8")

        def log(msg):
            line = f"{time.strftime('%H:%M:%S')} {msg}"
            print(line, flush=True)
            logf.write(line + "\n")
            logf.flush()
        if os.name == "nt":
            # keep Windows from idle-sleeping while this process runs (released on exit)
            import ctypes
            ctypes.windll.kernel32.SetThreadExecutionState(0x80000000 | 0x00000001)
        RUN_CONFIG.write_text(json.dumps({"b_full": args.b_full, "b_sub": args.b_sub}))
        log(f"run start: b_full={args.b_full} b_sub={args.b_sub} n_cpus={args.n_cpus}")
        run(args.datasets, args.b_full, args.b_sub, args.n_cpus, args.max_workers,
            args.r_worker_gb, args.reserve_gb, log)



# ---------------------------------------------------------------- summary (after the run)

FDR_CUT, LFC_CUT = 0.05, 1.0
PROTOCOL_ORDER = ["P1", "P2", "P3", "P4", "P5", "P6", "P7"]


def load_npz(ds: str, p: str, rep: str) -> tuple[np.ndarray, np.ndarray]:
    z = np.load(result_path(ds, p, rep))
    return dm.bh(z["pvalue"]), z["log2FC"].astype(float)


def collect_fp(datasets=None) -> tuple[pd.DataFrame, dict]:
    """Per replicate FP counts at both thresholds; also per-gene hit counts (FDR < 0.05)."""
    cfg = json.loads(RUN_CONFIG.read_text())
    rows, gene_hits = [], {}
    for ds, (_, protocols) in DATASETS.items():
        if datasets and ds not in datasets:
            continue
        for p in protocols:
            hits = None
            for rep in target_reps(ds, cfg["b_full"], cfg["b_sub"]):
                fdr, lfc = load_npz(ds, p, rep)
                sig = fdr < FDR_CUT
                rows.append({"dataset": ds, "protocol": p, "rep": rep,
                             "FP_fdr": int(sig.sum()),
                             "FP_fdr_lfc": int((sig & (np.abs(lfc) > LFC_CUT)).sum())})
                hits = sig.astype(int) if hits is None else hits + sig
            gene_hits[(ds, p)] = hits
    return pd.DataFrame(rows), gene_hits


def boot_ci(x: np.ndarray, n_boot: int = 2000, seed: int = SEED) -> tuple[float, float]:
    rng = np.random.default_rng(seed)
    x = np.asarray(x, float)
    means = rng.choice(x, size=(n_boot, len(x)), replace=True).mean(axis=1)
    return float(np.percentile(means, 2.5)), float(np.percentile(means, 97.5))


def real_counts(dataset_real: str | None, protocol: str) -> tuple[int, int]:
    """(n FDR < 0.05, n DEG) for real labels; None -> full-data Stage 3 tables."""
    if dataset_real is None:
        r = dm.load_protocol(protocol)
        return int((r["FDR"] < FDR_CUT).sum()), int(r["is_DEG"].sum())
    fdr, lfc = load_npz(dataset_real, protocol, "real")
    sig = fdr < FDR_CUT
    return int(sig.sum()), int((sig & (np.abs(lfc) > LFC_CUT)).sum())


def main_dataset(null: str, protocol: str) -> str:
    """Rank protocols P4/P6 are evaluated on the full sets, the rest on subsamples."""
    tag = "full" if protocol in RANK_SET else "sub"
    return f"{null}_{tag}"


def summarize_fp(fp: pd.DataFrame) -> pd.DataFrame:
    """One row per protocol x threshold: FP_A, FP_rand, FP_conf, FP_batch, FP_C, FDR_A."""
    rows = []
    for metric, thr in [("FP_fdr", "FDR<0.05"), ("FP_fdr_lfc", "FDR<0.05 & |log2FC|>1")]:
        for p in PROTOCOL_ORDER:
            row = {"protocol": p, "threshold": thr}
            for null in ["A", "Brand", "Bconf", "C"]:
                ds = main_dataset(null, p)
                x = fp[(fp["dataset"] == ds) & (fp["protocol"] == p)].sort_values("rep")[metric].to_numpy()
                lo, hi = boot_ci(x)
                key = {"A": "FP_A", "Brand": "FP_rand", "Bconf": "FP_conf", "C": "FP_C"}[null]
                row.update({key: x.mean(), f"{key}_median": float(np.median(x)), f"{key}_ci_lo": lo,
                            f"{key}_ci_hi": hi, f"{key}_frac_reps_any": float((x > 0).mean()),
                            f"{key}_dataset": ds, f"{key}_B": len(x)})
                if null == "Brand":
                    rand = x
                if null == "Bconf":
                    diff = x - rand  # same replicate: same samples, same number of fake normals
                    lo, hi = boot_ci(diff)
                    row.update({"FP_batch": diff.mean(), "FP_batch_ci_lo": lo, "FP_batch_ci_hi": hi})
            real_ds = None if p in RANK_SET else "real_subA"
            n_sig, n_deg = real_counts(real_ds, p)
            n_real = n_sig if metric == "FP_fdr" else n_deg
            row.update({"real_positives": n_real, "real_dataset": real_ds or "full (Stage 3)",
                        "FDR_A": row["FP_A"] / n_real if n_real else np.nan})
            rows.append(row)
    return pd.DataFrame(rows)


def batch_sensitive_genes(gene_hits: dict, genes: pd.DataFrame, index: pd.Index,
                          cut: float = 0.10) -> pd.DataFrame:
    """Per gene: share of B-conf (and B-rand) replicates in which it is significant
    (FDR < 0.05); batch-sensitive if the B-conf share is >= cut."""
    cfg = json.loads(RUN_CONFIG.read_text())
    out = pd.DataFrame(index=index)
    out.insert(0, "gene_name", genes.loc[index, "gene_name"].to_numpy())
    for p in PROTOCOL_ORDER:
        for null, tag in [("Bconf", "conf"), ("Brand", "rand")]:
            ds = main_dataset(null, p)
            b = cfg["b_full"] if ds.endswith("_full") else cfg["b_sub"]
            out[f"freq_{tag}_{p}"] = gene_hits[(ds, p)] / b
        out[f"batch_sensitive_{p}"] = out[f"freq_conf_{p}"] >= cut
        # stricter: frequent under confounded labels but not under random labels
        out[f"batch_specific_{p}"] = out[f"batch_sensitive_{p}"] & (out[f"freq_rand_{p}"] < cut)
    out["n_protocols_batch_sensitive"] = out[[f"batch_sensitive_{p}" for p in PROTOCOL_ORDER]].sum(axis=1)
    out.index.name = "gene_id"
    return out


def run_summary() -> dict:
    """Stage 4 tables and Fig 5 from the saved replicate results."""
    import pickle
    import plots
    fp, hits = collect_fp()
    fp.to_csv(TABLES / "perm_fp_per_replicate.csv", index=False)
    with open(PERM / "gene_hits.pkl", "wb") as fh:
        pickle.dump(hits, fh)
    summary = summarize_fp(fp)
    summary.to_csv(TABLES / "perm_summary.csv", index=False)
    idx = pd.read_parquet(PROCESSED / "counts_filtered.parquet", columns=[]).index
    genes = pd.read_parquet(PROCESSED / "genes.parquet")
    bs = batch_sensitive_genes(hits, genes, idx)
    bs.to_csv(TABLES / "batch_sensitive_genes.csv")
    diag = pd.read_csv(PERM / "label_diagnostics.csv")
    diag.groupby("null").agg(n=("n", "first"), fake_normal=("n_fake_normal", "mean"),
                             V_mean=("V_plate", "mean"), V_sd=("V_plate", "std"),
                             corr_with_truth_mean=("corr_with_truth", "mean"),
                             corr_with_truth_sd=("corr_with_truth", "std"))         .to_csv(TABLES / "perm_label_diagnostics.csv")
    fig = plots.fig5_false_positives(summary)
    return {"fp": fp, "summary": summary, "batch_sensitive": bs, "fig": fig}


if __name__ == "__main__":
    main()
