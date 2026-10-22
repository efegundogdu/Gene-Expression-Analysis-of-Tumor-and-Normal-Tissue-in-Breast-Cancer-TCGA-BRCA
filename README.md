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
- **ComBat-seq** runs in R 4.6.1 with Bioconductor sva 3.60.0 through Rscript
  (inmoose has no Windows wheel).
