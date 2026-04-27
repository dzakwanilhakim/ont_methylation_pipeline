#!/bin/bash
# ============================================================
# 03_modkit_pileup.sh
# Run modkit pileup on barcode01–barcode24 sorted BAMs
#
# Key flags:
#   --cpg             → CpG context only
#   --combine-strands → merge _F and _R alignment counts per CpG
#                       into one beta value (correct for buffered FASTA)
#   --only-tabs       → tab-delimited BED, no header comments
#   --ref             → required for CpG context detection
#
# Project structure:
#   myproject/
#   ├── data/
#   │   └── Buffered_5000_Fasta.fasta
#   ├── analysis/
#   │   ├── 03_aligned/    ← input
#   │   └── 04_pileup/     ← output
#   └── src/               ← this script lives here
# ============================================================
set -euo pipefail
#conda activate ont_methyl
trap 'echo "[ERROR] Script failed at line $LINENO: $BASH_COMMAND" >&2' ERR

# ── Project root ──────────────────────────────────────────────────────────────
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"

# ── Paths ─────────────────────────────────────────────────────────────────────
REF="$PROJECT_DIR/data/Buffered_5000_Fasta.fasta"

# DEMUX_DIR: same as used in 02_align.sh (only referenced for context; pileup reads from ALIGN_DIR)
DEMUX_DIR="$PROJECT_DIR/analysis/02_demuxed"
# Deep MinKNOW layout example (uncomment and edit if needed):
# DEMUX_DIR="$PROJECT_DIR/analysis/02_demuxed/<experiment>/<run_id>/bam_pass"

ALIGN_DIR="$PROJECT_DIR/analysis/03_aligned"
PILEUP_DIR="$PROJECT_DIR/analysis/04_pileup"

FIRST_BC=1
LAST_BC=24
THREADS=24

mkdir -p "$PILEUP_DIR"

echo "========================================"
echo " Step 3: modkit pileup"
echo " $(date)"
echo " Project   : $PROJECT_DIR"
echo " Reference : $REF"
echo " Barcodes  : $FIRST_BC – $LAST_BC"
echo " Threads   : $THREADS"
echo "========================================"
echo ""

for i in $(seq -w "$FIRST_BC" "$LAST_BC"); do
    BC_LABEL="barcode${i}"
    IN_BAM="$ALIGN_DIR/${BC_LABEL}.sorted.bam"
    OUT_BED="$PILEUP_DIR/${BC_LABEL}.bed"

    if [ ! -f "$IN_BAM" ]; then
        echo "[$(date)] WARNING: $IN_BAM not found — skipping $BC_LABEL"
        continue
    fi

    if [ "$(samtools view -c "$IN_BAM")" -eq 0 ]; then
        echo "[$(date)] WARNING: $BC_LABEL has 0 aligned reads — skipping"
        continue
    fi

    if [ -f "$OUT_BED" ] && [ "$(wc -l < "$OUT_BED")" -gt 0 ]; then
        echo "[$(date)] Skipping $BC_LABEL (valid pileup exists)"
        continue
    else
        rm -f "$OUT_BED"
    fi

    echo "[$(date)] modkit pileup: $BC_LABEL ..."

    modkit pileup \
        "$IN_BAM" \
        "$OUT_BED" \
        --ref "$REF" \
        --modified-bases C:m \
        --cpg \
        --combine-strands \
        --threads "$THREADS" \
        --log-filepath "$PILEUP_DIR/${BC_LABEL}_modkit.log"             

    LINES=$(wc -l < "$OUT_BED")
    echo "  → $BC_LABEL: $LINES CpG positions in pileup BED"
done

# ── Summary ───────────────────────────────────────────────────────────────────
echo ""
echo "=== Pileup Summary ==="
printf "  %-15s %10s\n" "Barcode" "CpG sites"
echo "  ─────────────────────────────"
for bed in "$PILEUP_DIR"/barcode*.bed; do
    [ -f "$bed" ] || continue
    bc=$(basename "$bed" .bed)
    lines=$(wc -l < "$bed")
    printf "  %-15s %10d\n" "$bc" "$lines"
done

echo ""
echo "[$(date)] Step 3 complete."
