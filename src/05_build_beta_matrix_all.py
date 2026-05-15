#!/usr/bin/env python3
"""
04_build_beta_matrix.py
=======================
Build the full CpG × barcode beta methylation matrix from modkit pileup BEDs,
covering every CpG along the entire 10002 bp buffered reference.

Design
------
The buffered FASTA contains two records per target CpG:
    cgXXXXXXXXX_F   →  5000 bp upstream + CG + 5000 bp downstream  (forward)
    cgXXXXXXXXX_R   →  reverse complement of _F                    (reverse)

Because _R is the reverse complement of _F, a CpG at position `p` on _F
corresponds to the same biological CpG at position `(10000 - p)` on _R.
We therefore mirror _R positions back to _F orientation and sum read
counts (Nmod_F + Nmod_R, Nvalid_F + Nvalid_R) for the same locus before
computing beta.

Row IDs encode position relative to the designed CpG center (pos 5000):
    cg10227331(0)       → the designed center CpG
    cg10227331(-760)    → CpG 760 bp upstream of center
    cg10227331(+2315)   → CpG 2315 bp downstream of center

Coverage filter is applied AFTER F+R combining, so a CpG with F=1 and R=5
will pass a threshold of 6.

Outputs (analysis/05_beta_matrix/):
    beta_matrix_full.csv       CpG × barcode, NaN = below threshold / no data
    beta_matrix_complete.csv   CpGs passing threshold in ALL barcodes
    beta_matrix_long.csv       tidy long format
    coverage_matrix.csv        combined Nvalid per CpG per barcode

modkit BED columns (--only-tabs --combine-strands at intra-record level):
  0:  chrom       e.g. cg10227331_F
  1:  start       0-based
  2:  end
  3:  mod_code    'm' = 5mC
  9:  Nvalid      M + U (combined-strands within this record)
  10: pct_mod
  11: Nmod
  12: Ncanon
"""

import os
import sys
import argparse
import pandas as pd
import numpy as np


REF_LEN = 10002          # buffered FASTA length: 5000 + CG + 5000
CENTER_POS = 5000        # 0-based BED start of the 'C' of the designed center CpG
MIRROR_BASE = 10000      # for _R → _F mirror: pos_F = MIRROR_BASE - pos_R


# ── Argument parsing ──────────────────────────────────────────────────────────
def parse_args():
    parser = argparse.ArgumentParser(
        description="Build full CpG × barcode beta matrix from modkit pileup BEDs"
    )
    parser.add_argument(
        "--project",
        default=None,
        help="Absolute path to project root. Defaults to one level up from this script."
    )
    parser.add_argument(
        "--min-cov",
        type=int,
        default=None,
        help="Minimum combined (F+R) Nvalid to keep a CpG. If omitted, prompt at runtime."
    )
    parser.add_argument(
        "--first-bc", type=int, default=1, help="First barcode number (default: 1)"
    )
    parser.add_argument(
        "--last-bc", type=int, default=24, help="Last barcode number (default: 24)"
    )
    return parser.parse_args()


# ── Load one barcode BED ──────────────────────────────────────────────────────
def load_pileup(bed_path: str, min_cov: int):
    """
    Load one modkit pileup BED, combine _F + _R counts at matching loci,
    apply coverage filter, and return per-CpG beta + coverage.

    Returns
    -------
    beta_series : pd.Series  index = CpG row label, values = beta (0–1)
    cov_series  : pd.Series  index = CpG row label, values = combined Nvalid
    """
    barcode = os.path.basename(bed_path).replace(".bed", "")

    try:
        df = pd.read_csv(bed_path, sep="\t", header=None, comment="#")
    except Exception as e:
        print(f"  ERROR reading {bed_path}: {e}")
        empty = pd.Series(dtype=float, name=barcode)
        return empty, empty

    if df.shape[1] < 13:
        print(f"  ERROR: unexpected BED format in {bed_path} ({df.shape[1]} cols)")
        empty = pd.Series(dtype=float, name=barcode)
        return empty, empty

    # Only keep the columns we need
    df = df.rename(columns={
        0: "chrom",
        1: "start",
        3: "mod_code",
        9: "Nvalid",
        11: "Nmod",
    })[["chrom", "start", "mod_code", "Nvalid", "Nmod"]]

    # 5mC CpG calls only
    df = df[df["mod_code"] == "m"].copy()
    df["start"] = pd.to_numeric(df["start"], errors="coerce").astype("Int64")
    df["Nvalid"] = pd.to_numeric(df["Nvalid"], errors="coerce").fillna(0).astype(int)
    df["Nmod"] = pd.to_numeric(df["Nmod"], errors="coerce").fillna(0).astype(int)
    df = df.dropna(subset=["start"])
    df["start"] = df["start"].astype(int)

    # Split strand suffix
    df["cpg"] = df["chrom"].str.replace(r"_[FR]$", "", regex=True)
    df["strand"] = df["chrom"].str.extract(r"_([FR])$", expand=False)

    # Map _R positions back to _F orientation:
    #   _F: keep pos as-is
    #   _R: pos_F_equivalent = MIRROR_BASE - pos_R
    df["pos_F"] = np.where(
        df["strand"] == "F",
        df["start"],
        MIRROR_BASE - df["start"],
    ).astype(int)

    # Drop any rows outside the valid window (shouldn't happen, but safe)
    df = df[(df["pos_F"] >= 0) & (df["pos_F"] < REF_LEN)]

    if df.empty:
        empty = pd.Series(dtype=float, name=barcode)
        return empty, empty

    # Sum F + R counts at the same (cpg, pos_F) locus
    combined = (
        df.groupby(["cpg", "pos_F"], as_index=False)
          .agg(Nmod=("Nmod", "sum"), Nvalid=("Nvalid", "sum"))
    )

    # Coverage filter on combined depth
    combined = combined[combined["Nvalid"] >= min_cov]

    if combined.empty:
        empty = pd.Series(dtype=float, name=barcode)
        return empty, empty

    # Build row labels: cg10227331(+/-offset)
    offset = combined["pos_F"] - CENTER_POS
    sign = np.where(offset > 0, "+", "")  # negatives already carry their minus sign
    combined["row_id"] = (
        combined["cpg"]
        + "("
        + sign + offset.astype(str)
        + ")"
    )

    beta = (combined["Nmod"] / combined["Nvalid"]).astype(float)
    beta.index = combined["row_id"].values
    beta.name = barcode

    cov = combined["Nvalid"].astype(float)
    cov.index = combined["row_id"].values
    cov.name = barcode

    return beta, cov


# ── Build matrix ──────────────────────────────────────────────────────────────
def build_matrix(pileup_dir, first_bc, last_bc, min_cov):
    bed_files = []
    missing = []

    for i in range(first_bc, last_bc + 1):
        bc = f"barcode{i:02d}"
        path = os.path.join(pileup_dir, f"{bc}.bed")
        if os.path.isfile(path):
            bed_files.append((bc, path))
        else:
            missing.append(bc)
            print(f"  WARNING: {path} not found")

    if not bed_files:
        print(f"ERROR: no BED files in {pileup_dir}", file=sys.stderr)
        sys.exit(1)

    print(f"\nFound {len(bed_files)}/{last_bc - first_bc + 1} barcode BEDs")
    if missing:
        print(f"Missing: {', '.join(missing)}")
    print(f"Window           : entire {REF_LEN} bp (all CpGs)")
    print(f"Center position  : {CENTER_POS} (0-based) → offset 0")
    print(f"Coverage filter  : combined F+R Nvalid ≥ {min_cov}")
    print()

    beta_list, cov_list = [], []

    for bc, path in bed_files:
        print(f"  Loading {bc} ...")
        beta_s, cov_s = load_pileup(path, min_cov)
        print(f"     → {len(beta_s)} CpGs passed coverage filter")
        beta_list.append(beta_s)
        cov_list.append(cov_s)

    # Union of all CpG row IDs, sorted by (cpg name, offset)
    all_ids = sorted(
        set().union(*[s.index for s in beta_list]),
        key=_sort_key,
    )

    beta_matrix = pd.concat(beta_list, axis=1).reindex(all_ids)
    cov_matrix = pd.concat(cov_list, axis=1).reindex(all_ids)
    beta_matrix.index.name = "CpG_ID"
    cov_matrix.index.name = "CpG_ID"

    return beta_matrix, cov_matrix


def _sort_key(row_id: str):
    """Sort row IDs by (cpg_name, signed_offset)."""
    cpg, _, rest = row_id.partition("(")
    offset = int(rest.rstrip(")"))
    return (cpg, offset)


# ── QC report ─────────────────────────────────────────────────────────────────
def qc_report(beta, cov, min_cov):
    n_cpg, n_bc = beta.shape

    print()
    print("=" * 60)
    print(" Beta Matrix QC Report")
    print("=" * 60)
    print(f"  Matrix shape          : {n_cpg} CpGs × {n_bc} barcodes")
    print(f"  CpGs with any data    : {beta.notna().any(axis=1).sum()}")
    print(f"  Fully observed CpGs   : {beta.notna().all(axis=1).sum()}")
    print(f"  CpGs ≥ {min_cov} cov in ALL bc : {(cov >= min_cov).fillna(False).all(axis=1).sum()}")
    print()

    # Missingness per barcode
    print("  Missing per barcode:")
    miss = beta.isna().sum().sort_values(ascending=False)
    for bc, n in miss.items():
        pct = 100 * n / n_cpg if n_cpg else 0
        flag = "  ← high missingness" if pct > 20 else ""
        print(f"    {bc:<15s}: {n:5d} ({pct:5.1f}%){flag}")

    # Coverage stats
    print()
    print("  Combined (F+R) coverage per barcode:")
    print(f"    {'Barcode':<15s} {'mean':>8} {'median':>8} {'min':>6} {'max':>6} {'n':>6}")
    print(f"    {'-'*15} {'-'*8} {'-'*8} {'-'*6} {'-'*6} {'-'*6}")
    for bc in cov.columns:
        v = cov[bc].dropna()
        if len(v) == 0:
            print(f"    {bc:<15s} {'NA':>8} {'NA':>8} {'NA':>6} {'NA':>6} {0:>6}")
            continue
        print(f"    {bc:<15s} {v.mean():8.1f} {v.median():8.0f} "
              f"{v.min():6.0f} {v.max():6.0f} {len(v):6d}")

    # Beta distribution
    flat = beta.to_numpy().flatten()
    flat = flat[~np.isnan(flat)]
    if flat.size == 0:
        print("\n  WARNING: no beta values for distribution analysis")
        print("=" * 60)
        return

    print()
    print(f"  Beta distribution (n={len(flat)}):")
    print(f"    mean   = {flat.mean():.3f}")
    print(f"    median = {np.median(flat):.3f}")
    print(f"    min    = {flat.min():.3f}")
    print(f"    max    = {flat.max():.3f}")

    bins = [0, 0.2, 0.4, 0.6, 0.8, 1.0001]
    labels = ["0–0.2", "0.2–0.4", "0.4–0.6", "0.6–0.8", "0.8–1.0"]
    counts, _ = np.histogram(flat, bins=bins)
    print("\n  Distribution:")
    for label, cnt in zip(labels, counts):
        pct = 100 * cnt / len(flat)
        bar = "█" * int(pct / 2)
        print(f"    {label:<8s}: {bar:<50s} {cnt:6d} ({pct:5.1f}%)")

    # Offset distribution (where do passing CpGs sit relative to center?)
    offsets = []
    for row_id in beta.index:
        if beta.loc[row_id].notna().any():
            offsets.append(_sort_key(row_id)[1])
    if offsets:
        offsets = np.array(offsets)
        print()
        print(f"  Position offsets (any-data CpGs, n={len(offsets)}):")
        print(f"    range: {offsets.min():+d} to {offsets.max():+d}")
        print(f"    at center (0)          : {(offsets == 0).sum()}")
        print(f"    upstream  (< 0)        : {(offsets < 0).sum()}")
        print(f"    downstream (> 0)       : {(offsets > 0).sum()}")
        print(f"    within ±3000 of center : {((offsets >= -3000) & (offsets <= 3000)).sum()}")

    print("=" * 60)


# ── Coverage threshold prompt ────────────────────────────────────────────────
def prompt_min_cov():
    while True:
        try:
            raw = input("Enter minimum coverage threshold (combined F+R Nvalid) [default 5]: ").strip()
            if raw == "":
                return 5
            val = int(raw)
            if val < 1:
                print("  Must be ≥ 1.")
                continue
            return val
        except ValueError:
            print("  Please enter an integer.")
        except EOFError:
            print("  No input available — using default 5.")
            return 5


# ── Main ──────────────────────────────────────────────────────────────────────
def main():
    args = parse_args()

    if args.project:
        project_dir = os.path.abspath(args.project)
    else:
        script_dir = os.path.dirname(os.path.abspath(__file__))
        project_dir = os.path.dirname(script_dir)

    pileup_dir = os.path.join(project_dir, "analysis", "04_pileup")
    out_dir = os.path.join(project_dir, "analysis", "05_beta_matrix")
    os.makedirs(out_dir, exist_ok=True)

    # Resolve coverage threshold
    min_cov = args.min_cov if args.min_cov is not None else prompt_min_cov()

    print("=" * 60)
    print(" ONT Methylation Beta Matrix Builder (full CpG, F+R combined)")
    print(f" {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 60)
    print(f"  Project dir : {project_dir}")
    print(f"  Pileup dir  : {pileup_dir}")
    print(f"  Output dir  : {out_dir}")
    print(f"  Barcodes    : {args.first_bc:02d} – {args.last_bc:02d}")
    print(f"  Min coverage: {min_cov}")

    beta, cov = build_matrix(pileup_dir, args.first_bc, args.last_bc, min_cov)

    if beta.empty:
        print("ERROR: beta matrix is empty. Check pileup inputs / coverage threshold.",
              file=sys.stderr)
        sys.exit(1)

    qc_report(beta, cov, min_cov)

    # Save outputs
    p_full = os.path.join(out_dir, "beta_matrix_full.csv")
    beta.to_csv(p_full, float_format="%.4f")
    print(f"\nSaved: {p_full}")

    complete_mask = beta.notna().all(axis=1)
    p_complete = os.path.join(out_dir, "beta_matrix_complete.csv")
    beta.loc[complete_mask].to_csv(p_complete, float_format="%.4f")
    print(f"Saved: {p_complete}  ({complete_mask.sum()} CpGs present in all barcodes)")

    p_cov = os.path.join(out_dir, "coverage_matrix.csv")
    cov.to_csv(p_cov, float_format="%.0f")
    print(f"Saved: {p_cov}")

    long_b = beta.reset_index().melt(id_vars="CpG_ID", var_name="barcode", value_name="beta")
    long_c = cov.reset_index().melt(id_vars="CpG_ID", var_name="barcode", value_name="Nvalid")
    long = (long_b.merge(long_c, on=["CpG_ID", "barcode"], how="left")
                  .dropna(subset=["beta"])
                  .sort_values(["CpG_ID", "barcode"])
                  .reset_index(drop=True))
    p_long = os.path.join(out_dir, "beta_matrix_long.csv")
    long.to_csv(p_long, index=False, float_format="%.4f")
    print(f"Saved: {p_long}  ({len(long)} rows)")

    print()
    print("Done.")


if __name__ == "__main__":
    main()
