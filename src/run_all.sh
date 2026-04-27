#!/bin/bash
# ============================================================
# run_all.sh
# Master script — runs the complete ONT methylation pipeline
#
# Project structure:
#   myproject/
#   ├── data/
#   │   ├── pod5/                     ← raw pod5 input
#   │   └── Buffered_5000_Fasta.fasta ← buffered reference
#   ├── analysis/                     ← all outputs land here
#   │   ├── 01_basecalled/
#   │   ├── 02_demuxed/
#   │   ├── 03_aligned/
#   │   ├── 04_pileup/
#   │   └── 05_beta_matrix/
#   └── src/                          ← all scripts live here
#       ├── run_all.sh                (this file)
#       ├── 00_setup_env.sh
#       ├── 01_basecall_demux.sh
#       ├── 02_align.sh
#       ├── 03_modkit_pileup.sh
#       └── 04_build_beta_matrix.py
#
# Usage:
#   bash src/run_all.sh               # run all steps
#   bash src/run_all.sh --from 2      # resume from step 2
#   bash src/run_all.sh --from 3 --min-cov 10
# ============================================================
set -euo pipefail

# ── Locate scripts ────────────────────────────────────────────────────────────
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"

# ── Defaults ──────────────────────────────────────────────────────────────────
FROM_STEP=1
MIN_COV=5

# ── Parse arguments ───────────────────────────────────────────────────────────
while [[ $# -gt 0 ]]; do
    case "$1" in
        --from)    FROM_STEP="$2"; shift 2 ;;
        --min-cov) MIN_COV="$2";   shift 2 ;;
        *)
            echo "Unknown argument: $1"
            echo "Usage: bash run_all.sh [--from STEP] [--min-cov N]"
            exit 1
            ;;
    esac
done

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }
hr()  { echo ""; echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"; }

hr
log " ONT Methylation Pipeline"
log " Project : $PROJECT_DIR"
log " Starting from step $FROM_STEP  |  min-cov=$MIN_COV"
hr

# ── Step 1: Basecall + Demux ──────────────────────────────────────────────────
if [ "$FROM_STEP" -le 1 ]; then
    hr
    log "STEP 1/4 — Basecalling + Demultiplexing (dorado, SQK-RBK114-96)"
    hr
    bash "$SCRIPT_DIR/01_basecall_demux.sh"
    log "✓ Step 1 done"
else
    log "Skipping Step 1 (--from $FROM_STEP)"
fi

# ── Step 2: Align ─────────────────────────────────────────────────────────────
if [ "$FROM_STEP" -le 2 ]; then
    hr
    log "STEP 2/4 — Alignment (minimap2 → barcode01–barcode24)"
    hr
    bash "$SCRIPT_DIR/02_align.sh"
    log "✓ Step 2 done"
else
    log "Skipping Step 2 (--from $FROM_STEP)"
fi

# ── Step 3: Pileup ────────────────────────────────────────────────────────────
if [ "$FROM_STEP" -le 3 ]; then
    hr
    log "STEP 3/4 — Methylation pileup (modkit)"
    hr
    bash "$SCRIPT_DIR/03_modkit_pileup.sh"
    log "✓ Step 3 done"
else
    log "Skipping Step 3 (--from $FROM_STEP)"
fi

# ── Step 4: Beta matrix ───────────────────────────────────────────────────────
if [ "$FROM_STEP" -le 4 ]; then
    hr
    log "STEP 4/4 — Build beta matrix (python)"
    hr
    conda activate ontmethyl_env
    python "$SCRIPT_DIR/04_build_beta_matrix.py" \
        --project "$PROJECT_DIR" \
        --min-cov "$MIN_COV"
    log "✓ Step 4 done"
else
    log "Skipping Step 4 (--from $FROM_STEP)"
fi

# ── Done ──────────────────────────────────────────────────────────────────────
hr
log "Pipeline complete."
echo ""
echo "  Final outputs:"
echo "  $PROJECT_DIR/analysis/05_beta_matrix/"
echo "    ├── beta_matrix_full.csv       (N_CpG × 24, NaN = low coverage)"
echo "    ├── beta_matrix_complete.csv   (CpGs present in all 24 barcodes)"
echo "    ├── beta_matrix_long.csv       (tidy format for ML / ggplot)"
echo "    └── coverage_matrix.csv        (Nvalid per CpG per barcode)"
hr
