"""Stage 1: load GDC STAR count files, map files to aliquot barcodes and build sample table."""

from __future__ import annotations

import glob
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd

SEED = 42

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
VERI = RAW / "veri"
PROCESSED = ROOT / "data" / "processed"

SAMPLE_TYPE_LABELS = {"01": "Primary Tumor", "11": "Solid Tissue Normal", "06": "Metastatic"}
CONDITION = {"01": "tumor", "11": "normal"}
STAR_SUMMARY_PREFIX = "N_"


# ---------------------------------------------------------------- counts

def find_count_files(veri_dir: Path = VERI) -> dict[str, Path]:
    """Return {file_id (folder UUID): path to the single STAR counts .tsv}."""
    files = {}
    for folder in sorted(p for p in veri_dir.iterdir() if p.is_dir()):
        tsvs = sorted(folder.glob("*.rna_seq.augmented_star_gene_counts.tsv"))
        if len(tsvs) != 1:
            raise ValueError(f"{folder.name}: expected 1 counts tsv, found {len(tsvs)}")
        files[folder.name] = tsvs[0]
    return files


def read_star_counts(path: Path) -> pd.DataFrame:
    """Read one GDC STAR augmented gene counts file, dropping the N_* summary rows."""
    df = pd.read_csv(path, sep="\t", comment="#",
                     usecols=["gene_id", "gene_name", "gene_type", "unstranded"])
    return df[~df["gene_id"].str.startswith(STAR_SUMMARY_PREFIX)].reset_index(drop=True)


def build_count_matrix(veri_dir: Path = VERI) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Build a genes x files integer matrix from the unstranded column.

    Returns (counts, genes): counts indexed by gene_id with file_id columns,
    genes with gene_id/gene_name/gene_type annotation.
    """
    files = find_count_files(veri_dir)
    genes = None
    columns = {}
    for file_id, path in files.items():
        df = read_star_counts(path)
        if genes is None:
            genes = df[["gene_id", "gene_name", "gene_type"]].copy()
        elif not df["gene_id"].equals(genes["gene_id"]):
            raise ValueError(f"{file_id}: gene_id order differs from the first file")
        columns[file_id] = df["unstranded"].to_numpy()
    counts = pd.DataFrame(columns, index=genes["gene_id"].to_numpy())
    counts.index.name = "gene_id"
    if (counts.values < 0).any() or not np.issubdtype(counts.values.dtype, np.integer):
        raise ValueError("counts must be non-negative integers")
    counts = counts.astype(np.int32)
    return counts, genes.set_index("gene_id")


# ---------------------------------------------------------------- metadata

def load_sample_sheet(raw_dir: Path = RAW) -> pd.DataFrame:
    paths = sorted(glob.glob(str(raw_dir / "gdc_sample_sheet*.tsv")))
    if len(paths) != 1:
        raise FileNotFoundError(f"expected one gdc_sample_sheet*.tsv, found {paths}")
    return pd.read_csv(paths[0], sep="\t")


def load_aliquots(raw_dir: Path = RAW) -> pd.DataFrame:
    """Map file_id -> aliquot barcode from metadata.repository*.json."""
    paths = sorted(glob.glob(str(raw_dir / "metadata.repository*.json")))
    if len(paths) != 1:
        raise FileNotFoundError(f"expected one metadata.repository*.json, found {paths}")
    with open(paths[0], encoding="utf-8") as fh:
        records = json.load(fh)
    rows = []
    for rec in records:
        ents = rec["associated_entities"]
        if len(ents) != 1 or ents[0]["entity_type"] != "aliquot":
            raise ValueError(f"{rec['file_id']}: unexpected associated_entities {ents}")
        rows.append({"file_id": rec["file_id"],
                     "aliquot_barcode": ents[0]["entity_submitter_id"]})
    return pd.DataFrame(rows)


def parse_barcode(barcode: str) -> dict:
    """Split a TCGA aliquot barcode, e.g. TCGA-A8-A09E-01A-11R-A00Z-07."""
    f = barcode.split("-")
    if len(f) != 7:
        raise ValueError(f"not an aliquot barcode: {barcode}")
    return {"tss": f[1], "patient": barcode[:12], "sample_type_code": f[3][:2],
            "vial": f[3][2:], "portion_analyte": f[4], "plate": f[5], "center": f[6]}


# ---------------------------------------------------------------- clinical

CLINICAL_COLS = {
    "demographic.age_at_index": "age",
    "diagnoses.ajcc_pathologic_stage": "stage",
    "demographic.race": "race",
    "demographic.vital_status": "vital_status",
}


def load_clinical(raw_dir: Path = RAW) -> tuple[pd.DataFrame, dict]:
    """Collapse clinical.tsv to one row per patient.

    For each field the first non-missing value is kept; the number of patients
    with conflicting non-missing values is returned for reporting.
    """
    path = raw_dir / "clinical" / "clinical.tsv"
    df = pd.read_csv(path, sep="\t", na_values=["'--"], low_memory=False)
    df = df[["cases.submitter_id", *CLINICAL_COLS]].rename(
        columns={"cases.submitter_id": "patient", **CLINICAL_COLS})
    grouped = df.groupby("patient", sort=True)
    conflicts = {c: int((grouped[c].nunique() > 1).sum()) for c in CLINICAL_COLS.values()}
    out = grouped.first()
    out["age"] = pd.to_numeric(out["age"], errors="coerce")
    return out.reset_index(), conflicts


# ---------------------------------------------------------------- samples

def build_samples(file_ids: list[str], raw_dir: Path = RAW) -> tuple[pd.DataFrame, dict]:
    """Build the filtered sample table plus a dict of QC numbers for the report."""
    report = {}
    sheet = load_sample_sheet(raw_dir)
    aliq = load_aliquots(raw_dir)

    # every folder must be in the sheet and the metadata json
    missing_sheet = sorted(set(file_ids) - set(sheet["File ID"]))
    missing_meta = sorted(set(file_ids) - set(aliq["file_id"]))
    if missing_sheet or missing_meta:
        raise ValueError(f"folders missing from sheet: {missing_sheet[:5]}, "
                         f"from metadata: {missing_meta[:5]}")
    report["n_files"] = len(file_ids)
    report["sheet_rows"] = len(sheet)

    df = aliq[aliq["file_id"].isin(file_ids)].merge(
        sheet[["File ID", "Sample ID", "Case ID", "Tissue Type", "Tumor Descriptor"]],
        left_on="file_id", right_on="File ID", how="left").drop(columns="File ID")
    df = pd.concat([df, pd.DataFrame([parse_barcode(b) for b in df["aliquot_barcode"]])], axis=1)

    # cross-checks between the two metadata sources
    report["sheet_sample_id_mismatch"] = int((df["aliquot_barcode"].str[:16] != df["Sample ID"]).sum())
    report["sheet_case_id_mismatch"] = int((df["patient"] != df["Case ID"]).sum())
    report["sample_type_counts_all"] = df["sample_type_code"].value_counts().sort_index().to_dict()
    report["centers"] = df["center"].value_counts().to_dict()

    # keep primary tumor and solid tissue normal only
    df = df[df["sample_type_code"].isin(CONDITION)].copy()
    df["condition"] = df["sample_type_code"].map(CONDITION)

    # one aliquot per patient and condition: alphabetically first barcode
    dup = df.groupby(["patient", "condition"])["aliquot_barcode"].transform("size") > 1
    report["patients_multi_tumor"] = int(df[dup & (df["condition"] == "tumor")]["patient"].nunique())
    report["patients_multi_normal"] = int(df[dup & (df["condition"] == "normal")]["patient"].nunique())
    report["dropped_duplicates"] = df[dup].sort_values("aliquot_barcode") \
        .groupby(["patient", "condition"]).apply(lambda g: list(g["aliquot_barcode"].iloc[1:])) \
        .explode().tolist() if dup.any() else []
    df = df.sort_values("aliquot_barcode").drop_duplicates(["patient", "condition"], keep="first")

    # paired subsets
    cond_per_patient = df.groupby("patient")["condition"].nunique()
    paired = set(cond_per_patient[cond_per_patient == 2].index)
    df["paired"] = df["patient"].isin(paired)
    plates = df[df["paired"]].pivot(index="patient", columns="condition", values="plate")
    same_plate = set(plates.index[plates["tumor"] == plates["normal"]])
    df["same_plate_paired"] = df["patient"].isin(same_plate)

    clinical, conflicts = load_clinical(raw_dir)
    report["clinical_conflicts"] = conflicts
    report["patients_without_clinical"] = int((~df["patient"].isin(clinical["patient"])).sum())
    df = df.merge(clinical, on="patient", how="left")

    cols = ["file_id", "aliquot_barcode", "patient", "condition", "sample_type_code", "tss",
            "plate", "center", "vial", "portion_analyte", "paired", "same_plate_paired",
            "age", "stage", "race", "vital_status"]
    df = df[cols].sort_values(["condition", "aliquot_barcode"],
                              ascending=[False, True]).reset_index(drop=True)
    return df, report


def summarize(samples: pd.DataFrame) -> dict:
    s = {}
    s["n_tumor"] = int((samples["condition"] == "tumor").sum())
    s["n_normal"] = int((samples["condition"] == "normal").sum())
    s["n_paired_patients"] = int(samples.loc[samples["paired"], "patient"].nunique())
    s["n_same_plate_paired_patients"] = int(samples.loc[samples["same_plate_paired"], "patient"].nunique())
    s["n_plates"] = int(samples["plate"].nunique())
    s["n_plates_with_normal"] = int(samples.loc[samples["condition"] == "normal", "plate"].nunique())
    s["n_tss"] = int(samples["tss"].nunique())
    tab = pd.crosstab(samples["plate"], samples["condition"])
    s["plates_normal_only"] = tab[(tab["tumor"] == 0)].to_dict("index")
    s["plates_tumor_only"] = int(((tab["normal"] == 0) & (tab["tumor"] > 0)).sum())
    return s


def main() -> None:
    PROCESSED.mkdir(parents=True, exist_ok=True)
    counts, genes = build_count_matrix()
    counts.to_parquet(PROCESSED / "counts.parquet")
    genes.to_parquet(PROCESSED / "genes.parquet")
    print(f"counts: {counts.shape[0]} genes x {counts.shape[1]} files")

    load_aliquots().to_csv(PROCESSED / "aliquots.csv", index=False)
    samples, report = build_samples(list(counts.columns))
    samples.to_csv(PROCESSED / "samples.csv", index=False)

    for k, v in {**report, **summarize(samples)}.items():
        print(f"{k}: {v}")


if __name__ == "__main__":
    main()
