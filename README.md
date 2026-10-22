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

**Start (or resume) the run** - the same command does both; finished work is skipped:

```bash
.venv/Scripts/python -I -u src/permutation.py run --b-full 100 --b-sub 50 --n-cpus 8 >> results/logs/perm_run.out 2>&1
```

Run it in the background (e.g. append `&` in bash, or a background task in Claude Code).
It writes progress to `results/logs/permutation.log` (one line per finished protocol run).
ComBat-seq runs in parallel R processes; their number is chosen at start so that total
RAM use stays below 70% of physical memory (`--max-workers` caps it, default 6).

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
