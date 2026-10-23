# Gene Expression Analysis of Tumor and Normal Tissue in Breast Cancer (TCGA-BRCA)

Work in progress. The full README (aim, data, protocols, how to run, findings,
limitations) is written in Stage 7. Method notes are collected here as decisions are made.

## Method notes

- **P7 winsorization (within condition).** P7 follows Yang et al. 2024 (Genome Biology
  25:282): counts are divided by DESeq2 size factors, each gene's normalized values above
  its 95th percentile are replaced by that percentile, then multiplied back by the size
  factors and rounded. The paper does not state whether the percentile is computed over
  all samples or per group. Computing it over all samples clipped the normals (~9% of
  samples) on down-regulated genes (ADH1B log2FC -3.57 vs -5.10 in P1), so the
  percentile is computed **within each condition**. The earlier global version is kept as
  `results/tables/de_P7_global.csv`. In permutation analyses the winsorization follows
  the permuted labels.
- **P5 permutations on the subsample.** ComBat-seq takes ~12 min per run on the full
  data (1208 samples), so P5 permutations, like the DESeq2 protocols, run on the
  balanced subsample instead of the full data. P3 and P5 share one ComBat-seq run per
  replicate.
- **Common log2FC sensitivity** (`results/tables/de_common_lfc/`). DEGs are re-called
  with each protocol's own FDR and a shared DESeq2 log2FC: P1 shrunken LFC for P1, P4,
  P7; P2 (~batch + condition) for P2, P6; P3 (ComBat-seq) for P3, P5.
- **Replicates.** Rank protocols on the full sets (P4, P6): B = 100. Subsample protocols
  (P1, P2, P3, P5, P7, plus P4/P6 on the subsamples as a check): B = 50, because of
  the DESeq2 and ComBat-seq cost (~6 h on 14 cores / 32 GB RAM).
- **Stage 4 null designs** (`src/permutation.py`, summary: `python src/permutation.py summary`).
  Null A shuffles tumor/normal labels over all samples. Null B uses tumors only (one per
  patient; every DEG is false): B-rand picks fake normals at random, B-conf mimics the real
  plate structure (~20% of fake normals on plates that become fake-normal-only, the rest on
  "mixed" plates at real normal fractions) and is redrawn until its plate x label Cramér's V
  is within 0.03 of the real data's V (0.552 full, 0.630 subsample). Both use the same
  samples and the same number of fake normals (9.35%); batch groups are rebuilt on the
  tumor-only sets with the Stage 2 rule. FP_batch = FP_conf - FP_rand, paired by replicate.
  Null C (secondary) shuffles labels within the 15 batch groups that hold both conditions.
  The original "restricted permutation over all batch groups" was dropped: single-condition
  groups keep their true labels, so permuted and true labels correlated at 0.29 (full) and
  0.44 (subsample) and real biology leaked into the null. Subsamples: Null A all 113
  normals + 339 random tumors; Null B 452 random tumors. FDR_A = mean FP_A / number of
  real-label positives on the same samples (subsample real run for P1-P3, P5, P7; full
  Stage 3 results for P4, P6). Error bars are 95% bootstrap CIs over replicates (seed 42).
- **Batch-sensitive genes** (`results/tables/batch_sensitive_genes.csv`): significant
  (FDR < 0.05) in >= 10% of B-conf replicates of a protocol. `batch_specific_*` additionally
  requires < 10% of B-rand replicates, separating batch-driven genes from genes that a
  method flags under any random labelling.
- **Limitations of the null designs.** Null B: tumor subtypes may be unevenly spread over
  plates, so part of FP_conf can come from biology tied to plates rather than technical batch.
  Null C: a residual label correlation remains (mean 0.07 full, 0.11 subsample) through
  differing tumor fractions between batch groups, and FP counts rise with it (Spearman
  0.83 for P1), so Null C mixes batch and biological signal.
- **P7 global vs within-condition winsorization** (`results/tables/p7_winsorization_comparison.csv`).
  A sensitivity variant "P7g" uses one global 95th percentile (otherwise the Yang et al.
  procedure). On the Null A subsample (B = 50, paired replicates) P7g cut false positives
  by 94% vs P1 at FDR < 0.05 (357 -> 21) and by 98.7% with |log2FC| > 1, in line with Yang
  et al. (98.2% fewer DEGs for DESeq2 at the 93rd percentile). Within-condition P7 cut them
  by only 3% (347). On real labels P7g finds fewer DEGs (3401 vs 3785 for P7 on the
  subsample). P7 in the main tables remains the within-condition version.
- **Stage 5 reference** (`src/benchmark.py`): the 89 same-plate tumor/normal pairs (178
  samples, 14 plates). Wilcoxon signed-rank on per-patient differences of log2(norm + 1)
  (log2FC = median difference) and PyDESeq2 `~ patient + condition`; reference = DEG
  (FDR < 0.05, |log2FC| > 1) in both, same direction. A protocol DEG is a true positive only
  with the same direction. Caveats: the 178 reference samples are also part of the full
  data the protocols use, and the Wilcoxon half of the reference uses a median-based
  log2FC, which favours the rank protocols' median log2FC; `benchmark_common_lfc.csv`
  repeats the benchmark with the shared DESeq2 log2FC.
- **UpSet plot** (Fig 6): intersection sizes come from `upsetplot.query`; the figure is
  drawn with matplotlib because upsetplot 0.9.0 plotting fails under pandas 3.
- **Stage 6 enrichment** (`src/downstream.py`): gseapy 1.3.1 `enrichr` in offline mode
  (hypergeometric test) on the Enrichr libraries GO_Biological_Process_2023 and
  KEGG_2021_Human, downloaded once to `data/processed/enrichr/`; background = the 15317
  tested genes. Top-20 terms = the 20 smallest adjusted p among terms with FDR < 0.05
  (KEGG up lists have only 8-12 such terms, so their Jaccard uses fewer terms).
- **Stage 6 classification**: nested stratified 5-fold CV (inner 5-fold grid search on
  ROC-AUC; L1-LR C in {0.01, 0.1, 1}; RF 300 trees, max_features in {sqrt, 0.02}),
  class_weight='balanced'. DEG selection is redone inside each outer training fold with
  the protocol's own test (P1, P2, P4, P6, P7). P3 and P5 are not used as in-fold
  selectors: ComBat-seq needs ~10 min per fold and there is no label-free way to adjust
  the held-out fold. Size factors use a median-of-ratios reference from the training fold
  only. ROC-AUC is at ceiling (>= 0.9994) for every protocol. When all C values tie on
  inner AUC, GridSearchCV keeps the first (C = 0.01), giving a very sparse model that
  misses normals at the 0.5 threshold in one fold (P7, P1_fixed: F1 normal ~0.72-0.75).
  `core_fixed` and `P1_fixed` use gene lists chosen on all samples before CV, so they are
  not nested estimates; they compare the two signatures under identical conditions.
- **Robust core** (`results/tables/robust_core.csv`): DEG in P2, P3, P5 and P6 with the
  same direction (2628 genes), minus genes significant in >= 5% of B-conf replicates of any
  of them (964 genes; 674 removed through P5, 346 only through P5) -> 1664 genes. The
  adipocyte markers ADH1B, FABP4 and LEP drop out because P5 flags them in 12-22% of B-conf
  replicates (0% for P6), i.e. the P5 (ComBat-seq + Wilcoxon) false positives drive most
  of the filtering. **Version 2** (`results/tables/robust_core_v2.csv`) keeps the
  four-protocol DEG/direction rule but applies the B-conf filter only to P2, P3 and P6:
  2010 genes (810 up, 1200 down; a superset of v1). ADH1B and FABP4 return; LEP stays out
  (6% of B-conf replicates in P2 and P3). The 346 added genes agree with the paired
  reference in 98.6% of cases (v1 overall: 94.9%).
- **ComBat-seq** runs in R 4.6.1 with Bioconductor sva 3.60.0 through Rscript
  (inmoose has no Windows wheel).

## Stopping and resuming the permutation run (Stage 4)

All commands run from the repository root. Python is the project venv
(`.venv/Scripts/python` on Windows, `.venv/bin/python` elsewhere). ComBat-seq needs
`Rscript` (R 4.6.1 with Bioconductor `sva`); it is found under `C:\Program Files\R`.

**Inputs that must already exist** (made by Stages 1-2 and the label step):
`data/processed/counts_filtered.parquet`, `data/processed/samples.csv` and the label sets
`data/processed/perm/labels_*.parquet`. If the label files are missing, recreate them
with `python src/permutation.py labels --B 100` (deterministic, seed 42; same labels every time).

**R packages.** Rscript is launched with `R_LIBS_USER` set to the project folder `.rlib/`
(gitignored) holding `sva` and its dependencies. This matters on Windows: R libraries
installed from inside a packaged app (such as the Claude desktop app) are redirected into
that app's private storage and are invisible to a scheduled task. If `.rlib/` is missing,
recreate it with `Rscript -e ".libPaths('.rlib'); install.packages('BiocManager', lib='.rlib',
repos='https://cloud.r-project.org'); BiocManager::install('sva', lib='.rlib', ask=FALSE, type='binary')"`.

**Recommended on Windows: start (or resume) as a scheduled task.** This keeps the run
alive when a terminal or Claude session closes (a run started from a Claude session was
killed when the session ended):

```powershell
powershell -NoProfile -ExecutionPolicy Bypass -File scripts\start_permutation_task.ps1
```

It registers the one-off task `TCGA_BRCA_permutation`, which runs
`scripts/run_permutation_task.cmd` hidden. The wrapper appends to
`results/logs/perm_task.log` (start time, runner exit code, final `status`) and deletes the
task when the runner exits - whether finished, stopped with the STOP file, or failed. Check
whether it is still running with `Get-ScheduledTask -TaskName TCGA_BRCA_permutation`
(no task = not running). If an R worker fails, the runner logs `ERROR: ComBat-seq failed`
in `results/logs/permutation.log` and stops instead of skipping P3/P5.

**Start (or resume) the run directly** - the same command does both; finished work is skipped:

```bash
.venv/Scripts/python -I -u src/permutation.py run --b-full 100 --b-sub 50 --n-cpus 8 >> results/logs/perm_run.out 2>&1
```

Run it in the background (e.g. append `&` in bash, or a background task in Claude Code).
It writes progress to `results/logs/permutation.log` (one line per finished protocol run).
ComBat-seq runs in parallel R processes; their number is chosen at start so that total
RAM use stays below 70% of physical memory (`--max-workers` caps it, default 6).
On Windows the runner asks the OS not to idle-sleep while it runs (released when it
exits); closing a laptop lid still sleeps the machine.

**Check progress** (safe at any time, also while running):

```bash
.venv/Scripts/python -I src/permutation.py status
```

It prints finished replicates per dataset x protocol and `total: done / target`, then
`complete` or `incomplete`. Target is 2216 protocol runs (8 full-data rank cells x 100,
28 subsample cells x 50, plus 16 real-label runs).

**Stop safely:** create the file `data/processed/perm/STOP`:

```bash
touch data/processed/perm/STOP
```

The runner finishes the protocol run it is on (at most ~1 min), terminates its R
workers, logs `STOP file found ... exiting` and exits. The STOP file is removed
automatically by the next `run`. If the process was killed hard instead (crash, reboot,
closed terminal), nothing is corrupted - results are written atomically - but leftover R
processes may keep running; stop them first (PowerShell: `Get-Process Rscript | Stop-Process`).

**Resume:** run the start command again. Each protocol result is one file,
`data/processed/perm/results/<dataset>/<protocol>/<rep>.npz` (per-gene p-values and
log2FC), saved right after it finishes; only missing files are computed. ComBat-seq output
that was not yet used when the run stopped is recomputed.

**When `status` says `complete`:** continue with the Stage 4 summary (FP counts, FDR
estimates, batch-sensitive genes, Fig 5).
