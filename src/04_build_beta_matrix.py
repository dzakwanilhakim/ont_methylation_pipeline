#!/usr/bin/env python3
"""
04_build_beta_matrix.py
=======================
Build the final beta methylation matrix: N_CpGs × 24 barcodes

Filters modkit pileup BED to position 5001 only (0-based: 5000)
— the true CpG center in your 5000-CG-5000 = 10002 bp buffered FASTA.

Project structure expected:
    myproject/
    ├── data/
    │   └── Buffered_5000_Fasta.fasta
    ├── analysis/
    │   ├── 04_pileup/          ← input BED files
    │   └── 05_beta_matrix/     ← output CSVs (created automatically)
    └── src/
        └── 04_build_beta_matrix.py  ← this script

Run from anywhere:
    python src/04_build_beta_matrix.py
    python src/04_build_beta_matrix.py --project /absolute/path/to/myproject
    python src/04_build_beta_matrix.py --min-cov 10

modkit BED columns (--only-tabs --combine-strands):
  Col 0:  chrom      e.g. CG00001_F  (strip _F/_R → CpG ID)
  Col 1:  start      0-based position
  Col 2:  end
  Col 3:  mod_code   'm' = 5mC
  Col 4:  score
  Col 5:  strand     '.' when --combine-strands used
  Col 6:  start (dup)
  Col 7:  end (dup)
  Col 8:  color
  Col 9:  Nvalid     total valid reads (M + U)
  Col 10: pct_mod    0–100  (= beta × 100)
  Col 11: Nmod       methylated read count
  Col 12: Ncanon     unmethylated read count
  Col 13: Nother
  Col 14: Ndel
  Col 15: Nfail
  Col 16: Ndiff
  Col 17: Nnocall

Outputs:
    05_beta_matrix/beta_matrix_full.csv       N_CpG × 24 (NaN = low coverage)
    05_beta_matrix/beta_matrix_complete.csv   CpGs with data in ALL 24 barcodes
    05_beta_matrix/beta_matrix_long.csv       tidy long-format for ML / plotting
    05_beta_matrix/coverage_matrix.csv        Nvalid per CpG per barcode (QC)
"""

import os
import sys
import glob
import argparse
import pandas as pd
import numpy as np

# ── Argument parsing ──────────────────────────────────────────────────────────
def parse_args():
    parser = argparse.ArgumentParser(
        description="Build beta methylation matrix from modkit pileup BEDs"
    )
    parser.add_argument(
        "--project",
        default=None,
        help=(
            "Absolute path to project root (myproject/). "
            "Defaults to two levels up from this script's location."
        )
    )
    parser.add_argument(
        "--min-cov",
        type=int,
        default=5,
        help="Minimum valid read coverage to include a CpG call (default: 5)"
    )
    parser.add_argument(
        "--target-pos",
        type=int,
        default=5000,
        help=(
            "0-based BED start position of the true CpG in the buffered FASTA. "
            "Buffer=5000bp → default 5000."
        )
    )
    parser.add_argument(
        "--first-bc", type=int, default=1,  help="First barcode number (default: 1)"
    )
    parser.add_argument(
        "--last-bc",  type=int, default=24, help="Last barcode number (default: 24)"
    )
    return parser.parse_args()


# ── BED column names ──────────────────────────────────────────────────────────
MODKIT_COLS = [
    "chrom", "start", "end", "mod_code", "score", "strand",
    "start2", "end2", "color", "Nvalid", "pct_mod",
    "Nmod", "Ncanon", "Nother", "Ndel", "Nfail", "Ndiff", "Nnocall"
]


# ── Helper: strip _F / _R suffix ─────────────────────────────────────────────
# ── Helper: strip strand suffix (_F / _R) ────────────────────────────────────
def cpg_id(chrom: str) -> str:
    """
    CG00001_F  →  CG00001
    CG00002_R  →  CG00002
    CG00003_FWD → CG00003 (safe handling)
    CG00004     → CG00004
    """
    if "_" in chrom:
        return chrom.rsplit("_", 1)[0]
    return chrom


# ── Load one barcode BED ──────────────────────────────────────────────────────
def load_pileup(bed_path: str, target_pos: int, min_cov: int):
    """
    Returns:
        beta_series : pd.Series  index=CpG_ID, values=beta (0–1)
        cov_series  : pd.Series  index=CpG_ID, values=Nvalid
    """
    barcode = os.path.basename(bed_path).replace(".bed", "")

    try:
        df = pd.read_csv(
            bed_path,
            sep="\t",
            header=None,
            comment="#"
        )
    except Exception as e:
        print(f"  ERROR reading {bed_path}: {e}")
        return pd.Series(dtype=float, name=barcode), pd.Series(dtype=float, name=barcode)

    #  ADD HERE
    if df.shape[1] < 13:
        print(f"ERROR: Unexpected BED format in {bed_path}")
        return pd.Series(dtype=float, name=barcode), pd.Series(dtype=float, name=barcode)

    # THEN continue
    df = df.rename(columns={
        0: "chrom",
        1: "start",
        3: "mod_code",
        9: "Nvalid",
        11: "Nmod",
        12: "Ncanon",
        10: "pct_mod"
    })
    # Ensure numeric
    df["start"] = pd.to_numeric(df["start"], errors="coerce")
    df["Nvalid"] = pd.to_numeric(df["Nvalid"], errors="coerce")
    df["Nmod"] = pd.to_numeric(df["Nmod"], errors="coerce")

    # Keep CpG only
    df = df[df["mod_code"] == "m"]

    # Keep ONLY target CpG position (your design is single-point)
    #df = df.loc[(df["start"] - target_pos).abs() <= 2]

    # Coverage filter
    df = df[df["Nvalid"] >= min_cov]

    if df.empty:
        return pd.Series(dtype=float, name=barcode), pd.Series(dtype=float, name=barcode)

    # CpG ID cleanup
    df["cpg_id"] = df["chrom"].str.replace(r"_[FR]$", "", regex=True)

    # IMPORTANT FIX:
    # do NOT .first() → you must SUM counts
    grouped = (
        df.groupby("cpg_id")
        .apply(lambda x: x.loc[x["Nvalid"].idxmax()])
    )
    
    beta_series = (grouped["Nmod"] / grouped["Nvalid"]).rename(barcode)
    cov_series  = grouped["Nvalid"].rename(barcode)

    return beta_series, cov_series


# ── Build matrix ──────────────────────────────────────────────────────────────
def build_matrix(pileup_dir: str, first_bc: int, last_bc: int,
                 target_pos: int, min_cov: int):

    # ── Collect BED files safely ──────────────────────────────────────────────
    bed_files = []
    missing_barcodes = []

    for i in range(first_bc, last_bc + 1):
        bc_label = f"barcode{i:02d}"
        bed_path = os.path.join(pileup_dir, f"{bc_label}.bed")

        if os.path.isfile(bed_path):
            bed_files.append((bc_label, bed_path))
        else:
            missing_barcodes.append(bc_label)
            print(f"  WARNING: {bed_path} not found — missing from matrix")

    if not bed_files:
        print(f"ERROR: No BED files found in {pileup_dir}", file=sys.stderr)
        sys.exit(1)

    print(f"\nFound {len(bed_files)}/{last_bc - first_bc + 1} barcode BED files")
    if missing_barcodes:
        print(f"Missing barcodes: {', '.join(missing_barcodes)}")

    print(f"Target position : {target_pos}")
    print(f"Minimum coverage: {min_cov} reads\n")

    beta_list = []
    cov_list  = []
    all_cpgs  = set()

    # ── Load each barcode ─────────────────────────────────────────────────────
    for bc_label, bed in bed_files:
        print(f"  Loading {bc_label} ...")
        beta_s, cov_s = load_pileup(bed, target_pos, min_cov)

        beta_s.name = bc_label
        cov_s.name  = bc_label

        beta_list.append(beta_s)
        cov_list.append(cov_s)

        all_cpgs.update(beta_s.index)

    # ── FIX: enforce consistent CpG ordering across all samples ───────────────
    all_cpgs = sorted(all_cpgs)

    beta_matrix = pd.concat(beta_list, axis=1).reindex(all_cpgs)
    cov_matrix  = pd.concat(cov_list, axis=1).reindex(all_cpgs)

    beta_matrix.index.name = "CpG_ID"
    cov_matrix.index.name  = "CpG_ID"

    return beta_matrix, cov_matrix


# ── QC report ─────────────────────────────────────────────────────────────────
def qc_report(beta_matrix: pd.DataFrame, cov_matrix: pd.DataFrame,
              min_cov: int) -> None:

    n_cpg, n_bc = beta_matrix.shape

    # ── CORE STATS ────────────────────────────────────────────────────────────
    any_data_mask = beta_matrix.notna().any(axis=1)
    full_data_mask = beta_matrix.notna().all(axis=1)

    n_any     = any_data_mask.sum()
    n_complete = full_data_mask.sum()

    print()
    print("=" * 60)
    print(" Beta Matrix QC Report")
    print("=" * 60)

    print(f"  Matrix shape          : {n_cpg} CpGs × {n_bc} barcodes")
    print(f"  CpGs with any data    : {n_any}")
    print(f"  Fully observed CpGs   : {n_complete}")

    # ✔ FIX: real biological meaning of completeness
    min_cov_mask = (cov_matrix >= min_cov).fillna(False)
    high_cov_cpgs = min_cov_mask.all(axis=1).sum()

    print(f"  CpGs ≥ {min_cov} cov   : {high_cov_cpgs} (all barcodes)")
    print()

    # ── MISSING DATA ──────────────────────────────────────────────────────────
    print("  Missing values per barcode:")
    missing = beta_matrix.isna().sum().sort_values(ascending=False)

    for bc, n in missing.items():
        pct = 100 * n / n_cpg if n_cpg > 0 else 0
        flag = "  ← high missingness" if pct > 20 else ""
        print(f"    {bc:<15s}: {n:5d} ({pct:5.1f}%){flag}")

    # ── COVERAGE STATS ────────────────────────────────────────────────────────
    print()
    print("  Coverage per barcode (Nvalid):")
    print(f"    {'Barcode':<15s} {'mean':>8} {'median':>8} {'min':>6} {'max':>6} {'n':>6}")
    print(f"    {'-'*15} {'-'*8} {'-'*8} {'-'*6} {'-'*6} {'-'*6}")

    for bc in cov_matrix.columns:
        vals = cov_matrix[bc].dropna()

        if len(vals) == 0:
            print(f"    {bc:<15s} {'NA':>8} {'NA':>8} {'NA':>6} {'NA':>6} {0:>6}")
            continue

        print(
            f"    {bc:<15s} "
            f"{vals.mean():8.1f} "
            f"{vals.median():8.0f} "
            f"{vals.min():6.0f} "
            f"{vals.max():6.0f} "
            f"{len(vals):6d}"
        )

    # ── BETA DISTRIBUTION ─────────────────────────────────────────────────────
    flat = beta_matrix.to_numpy().flatten()
    flat = flat[~np.isnan(flat)]

    if flat.size == 0:
        print("\n  WARNING: No beta values available for distribution analysis")
        print("=" * 60)
        return

    print()
    print(f"  Beta distribution (n={len(flat)}):")
    print(f"    mean   = {flat.mean():.3f}")
    print(f"    median = {np.median(flat):.3f}")
    print(f"    min    = {flat.min():.3f}")
    print(f"    max    = {flat.max():.3f}")

    bins   = [0, 0.2, 0.4, 0.6, 0.8, 1.0001]
    labels = ["0–0.2", "0.2–0.4", "0.4–0.6", "0.6–0.8", "0.8–1.0"]

    counts, _ = np.histogram(flat, bins=bins)

    print("\n  Distribution:")
    for label, cnt in zip(labels, counts):
        pct = 100 * cnt / len(flat)
        bar = "█" * int(pct / 2)
        print(f"    {label:<8s}: {bar:<50s} {cnt:6d} ({pct:5.1f}%)")

    print("=" * 60)


# ── Main ──────────────────────────────────────────────────────────────────────
def main():
    args = parse_args()

    # ── Resolve project root ──────────────────────────────────────────────────
    if args.project:
        project_dir = os.path.abspath(args.project)
    else:
        script_dir  = os.path.dirname(os.path.abspath(__file__))
        project_dir = os.path.dirname(script_dir)

    pileup_dir = os.path.join(project_dir, "analysis", "04_pileup")
    out_dir    = os.path.join(project_dir, "analysis", "05_beta_matrix")
    os.makedirs(out_dir, exist_ok=True)

    # ── Header ────────────────────────────────────────────────────────────────
    print("=" * 60)
    print(" ONT Methylation Beta Matrix Builder")
    print(f" {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 60)
    print(f"  Project dir : {project_dir}")
    print(f"  Pileup dir  : {pileup_dir}")
    print(f"  Output dir  : {out_dir}")
    print(f"  Barcodes    : {args.first_bc:02d} – {args.last_bc:02d}")
    print()

    # ── Build matrix ──────────────────────────────────────────────────────────
    beta_matrix, cov_matrix = build_matrix(
        pileup_dir=pileup_dir,
        first_bc=args.first_bc,
        last_bc=args.last_bc,
        target_pos=args.target_pos,
        min_cov=args.min_cov,
    )

    # ── SAFETY CHECK ─────────────────────────────────────────────────────────
    if beta_matrix.empty:
        print("ERROR: Beta matrix is empty. Check pileup input.", file=sys.stderr)
        sys.exit(1)

    # ── QC ────────────────────────────────────────────────────────────────────
    qc_report(beta_matrix, cov_matrix, args.min_cov)

    # ── Save outputs ─────────────────────────────────────────────────────────

    # Full matrix
    out_full = os.path.join(out_dir, "beta_matrix_full.csv")
    beta_matrix.to_csv(out_full, float_format="%.4f")
    print(f"\nSaved: {out_full}")

    # COMPLETE CpGs = strictly all barcodes non-NaN
    complete_mask = beta_matrix.notna().all(axis=1)
    beta_complete = beta_matrix.loc[complete_mask]

    out_complete = os.path.join(out_dir, "beta_matrix_complete.csv")
    beta_complete.to_csv(out_complete, float_format="%.4f")
    print(f"Saved: {out_complete}  ({len(beta_complete)} CpGs fully complete)")

    # Coverage matrix
    out_cov = os.path.join(out_dir, "coverage_matrix.csv")
    cov_matrix.to_csv(out_cov, float_format="%.0f")
    print(f"Saved: {out_cov}")

    # ── SAFE long format conversion ──────────────────────────────────────────
    beta_long = beta_matrix.reset_index().melt(
        id_vars="CpG_ID",
        var_name="barcode",
        value_name="beta"
    )

    cov_long = cov_matrix.reset_index().melt(
        id_vars="CpG_ID",
        var_name="barcode",
        value_name="Nvalid"
    )

    # SAFE merge (alignment guaranteed by CpG_ID + barcode)
    beta_long = (
        beta_long
        .merge(cov_long, on=["CpG_ID", "barcode"], how="left")
        .dropna(subset=["beta"])
        .sort_values(["CpG_ID", "barcode"])
        .reset_index(drop=True)
    )

    out_long = os.path.join(out_dir, "beta_matrix_long.csv")
    beta_long.to_csv(out_long, index=False, float_format="%.4f")
    print(f"Saved: {out_long}  ({len(beta_long)} rows)")

    # ── Summary ──────────────────────────────────────────────────────────────
    print()
    print("Done. Outputs:")
    print(f"  {out_dir}/")
    print("    ├── beta_matrix_full.csv")
    print("    ├── beta_matrix_complete.csv")
    print("    ├── beta_matrix_long.csv")
    print("    └── coverage_matrix.csv")

if __name__ == "__main__":
    main()
