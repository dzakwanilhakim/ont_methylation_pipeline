#!/bin/bash
# ============================================================
# 02_align.sh
# Align barcode01–barcode24 BAMs to Buffered_5000_Fasta.fasta
#
# CRITICAL minimap2 flags:
#   -y              → carry MM/ML methylation tags into output BAM
#   --secondary=no  → one best hit per read; prevents MAPQ=0 caused
#                     by duplicate _F and _R targets in buffered FASTA
#   -ax map-ont     → Oxford Nanopore R10.4.1 long-read preset
#
# Project structure:
#   myproject/
#   ├── data/
#   │   ├── pod5/
#   │   └── Buffered_5000_Fasta.fasta
#   ├── analysis/
#   │   ├── 02_demuxed/    ← input
#   │   └── 03_aligned/    ← output
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

# DEMUX_DIR: set to the bam_pass folder produced by dorado demux.
# If your MinKNOW/dorado output uses a deep subfolder structure, override here.
# Default: flat layout produced by 01_basecall_demux.sh
DEMUX_DIR="$PROJECT_DIR/analysis/02_demuxed"
# Deep MinKNOW layout example (uncomment and edit if needed):
# DEMUX_DIR="$PROJECT_DIR/analysis/02_demuxed/<experiment>/<run_id>/bam_pass"

ALIGN_DIR="$PROJECT_DIR/analysis/03_aligned"

# ── Kit & barcodes ────────────────────────────────────────────────────────────
KIT="SQK-RBK114-96"
FIRST_BC=1
LAST_BC=24
THREADS=24

mkdir -p "$ALIGN_DIR"

echo "========================================"
echo " Step 2: Alignment"
echo " $(date)"
echo " Project   : $PROJECT_DIR"
echo " Reference : $REF"
echo " Barcodes  : $FIRST_BC – $LAST_BC"
echo " Threads   : $THREADS"
echo "========================================"

# ── Index reference ───────────────────────────────────────────────────────────
if [ ! -f "${REF}.fai" ]; then
    echo "[$(date)] Indexing reference..."
    samtools faidx "$REF" 2>/dev/null
fi

N_CONTIGS=$(grep -c '>' "$REF")
echo "[$(date)] Reference contigs: $N_CONTIGS  (expect 2 × N_CpGs: _F + _R per target)"
echo ""

# ── Align barcode01–barcode24 ─────────────────────────────────────────────────
for i in $(seq -w "$FIRST_BC" "$LAST_BC"); do
    BC_LABEL="barcode${i}"

    # Find BAM automatically
    IN_BAM=$(ls "$DEMUX_DIR/${BC_LABEL}"/*.bam 2>/dev/null | head -n 1)

    if [ -z "$IN_BAM" ]; then
        echo "[$(date)] WARNING: No BAM found for $BC_LABEL — skipping"
        continue
    fi

    OUT_BAM="$ALIGN_DIR/${BC_LABEL}.sorted.bam"

    # Skip if valid output exists
    if [ -f "$OUT_BAM" ] && [ "$(samtools view -c "$OUT_BAM" 2>/dev/null)" -gt 0 ]; then
        echo "[$(date)] Skipping $BC_LABEL (valid output exists)"
        continue
    else
        rm -f "$OUT_BAM" "$OUT_BAM.bai"
    fi

    # Skip empty input BAM
    if [ "$(samtools view -c "$IN_BAM" 2>/dev/null)" -eq 0 ]; then
        echo "[$(date)] WARNING: $BC_LABEL has 0 reads — skipping"
        continue
    fi

    echo "[$(date)] Aligning $BC_LABEL ..."

    # FIX: use dorado aligner (preserves MM/ML, no SAM/FASTQ issues)
    dorado aligner \
        -t "$THREADS" \
        "$REF" \
        "$IN_BAM" \
    | samtools sort -@ "$THREADS" -o "$OUT_BAM" -

    samtools index "$OUT_BAM"

    TOTAL=$(samtools view -c "$OUT_BAM" 2>/dev/null)
    MAPPED=$(samtools view -c -F 4 "$OUT_BAM" 2>/dev/null)
    UNMAPPED=$(samtools view -c -f 4 "$OUT_BAM" 2>/dev/null)

    if [ "$TOTAL" -gt 0 ]; then
        PCT=$(awk "BEGIN {printf \"%.1f\", 100*$MAPPED/$TOTAL}")
    else
        PCT="0.0"
    fi

    echo "  → $BC_LABEL: $MAPPED/$TOTAL mapped ($PCT%), $UNMAPPED unmapped"
done

# ── Alignment summary ─────────────────────────────────────────────────────────
echo ""
echo "=== Alignment Summary ==="
printf "%-15s %10s %10s %10s\n" "Barcode" "Total" "Mapped" "Pct"
echo "  ─────────────────────────────────────────────"

for bam in "$ALIGN_DIR"/barcode*.sorted.bam; do
    [ -f "$bam" ] || continue

    total=$(samtools view -c "$bam" 2>/dev/null)
    [ "$total" -eq 0 ] && continue   # 🔥 skip empty BAMs

    bc=$(basename "$bam" .sorted.bam)
    mapped=$(samtools view -c -F 4 "$bam" 2>/dev/null)

    pct=$(awk "BEGIN {printf \"%.1f%%\", 100*$mapped/$total}")

    printf "  %-13s %10d %10d %10s\n" "$bc" "$total" "$mapped" "$pct"
done

echo ""
echo "[$(date)] Step 2 complete."
