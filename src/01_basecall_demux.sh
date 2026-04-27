#!/bin/bash
# ============================================================
# 01_basecall_demux.sh
# Basecall all pod5 files with Dorado (SUP + 5mCG_5hmCG model)
# then demultiplex into per-barcode BAMs
#
# Kit    : SQK-RBK114-96
# Barcodes used: barcode01 – barcode24
#
# Project structure:
#   myproject/
#   ├── data/
#   │   ├── pod5/          ← raw pod5 input
#   │   └── Buffered_5000_Fasta.fasta
#   ├── analysis/
#   │   ├── 01_basecalled/
#   │   ├── 02_demuxed/
#   │   ├── 03_aligned/
#   │   ├── 04_pileup/
#   │   └── 05_beta_matrix/
#   └── src/               ← this script lives here
# ============================================================
set -euo pipefail
# Trap any error and print the line number + command that failed
trap 'echo "[ERROR] Script failed at line $LINENO: $BASH_COMMAND" >&2' ERR
#conda activate ont_methyl

# ── Project root (one level up from src/) ────────────────────────────────────
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"

# ── Paths ─────────────────────────────────────────────────────────────────────
POD5_DIR="$PROJECT_DIR/data/pod5"
BASECALL_DIR="$PROJECT_DIR/analysis/01_basecalled"
DEMUX_DIR="$PROJECT_DIR/analysis/02_demuxed"
# Deep MinKNOW layout example (uncomment and edit if needed):
# DEMUX_DIR="$PROJECT_DIR/analysis/02_demuxed/<experiment>/<run_id>/bam_pass"

# ── Dorado binary — UPDATE THIS PATH ─────────────────────────────────────────
DORADO="dorado"           # e.g. ~/tools/dorado-0.9.1-linux-x64/bin/dorado
DORADO_MODEL="sup,5mCG_5hmCG"     # SUP basecalling + 5mC/5hmC methylation model

# ── Kit & barcodes ────────────────────────────────────────────────────────────
KIT="SQK-RBK114-96"               # 96-well kit; demux filters to barcode01-24 below
FIRST_BC=1
LAST_BC=24

# ── Setup ─────────────────────────────────────────────────────────────────────
mkdir -p "$BASECALL_DIR" "$DEMUX_DIR"

echo "========================================"
echo " Step 1: Basecalling + Demultiplexing"
echo " $(date)"
echo " Project  : $PROJECT_DIR"
echo " pod5 dir : $POD5_DIR"
echo " Kit      : $KIT  (barcodes $FIRST_BC–$LAST_BC)"
echo "========================================"

echo "[$(date)] Counting pod5 files..."
N_POD5=$(find "$POD5_DIR" -name "*.pod5" | wc -l)
echo "  Found $N_POD5 pod5 files"

# ── Basecalling ───────────────────────────────────────────────────────────────
echo ""
echo "[$(date)] === BASECALLING ==="

$DORADO basecaller \
    "$DORADO_MODEL" \
    "$POD5_DIR" \
    --kit-name "$KIT" \
    --recursive \
    --no-trim \
    --min-qscore 8 \
    2> "$BASECALL_DIR/basecall.log" \
    > "$BASECALL_DIR/calls.bam"

echo "[$(date)] Basecalling done."
echo "  Output  : $BASECALL_DIR/calls.bam"
echo "  Reads   : $(samtools view -c $BASECALL_DIR/calls.bam)"

# ── Demultiplexing ────────────────────────────────────────────────────────────
echo ""
echo "[$(date)] === DEMULTIPLEXING ==="

$DORADO demux \
    "$BASECALL_DIR/calls.bam" \
    --kit-name "$KIT" \
    --output-dir "$DEMUX_DIR" \
    --no-trim \
    --threads 24 \
    2> "$DEMUX_DIR/demux.log"

echo "[$(date)] Demux done."

# ── Report: only barcodes 01–24 ───────────────────────────────────────────────


echo "=== Read counts: barcode01–barcode24 ==="
TOTAL_USED=0
for i in $(seq -w "$FIRST_BC" "$LAST_BC"); do
    # Glob: works regardless of run ID in filename
    BAM=$(ls "$DEMUX_DIR/barcode${i}/"*.bam 2>/dev/null | head -n 1)
    if [ -n "$BAM" ] && [ -f "$BAM" ]; then
        COUNT=$(samtools view -c "$BAM")
        TOTAL_USED=$((TOTAL_USED + COUNT))
        printf "  barcode%s : %d reads\n" "$i" "$COUNT"
    else
        printf "  barcode%s : MISSING\n" "$i"
    fi
done

UNCLASS_BAM=$(ls "$DEMUX_DIR/unclassified/"*.bam 2>/dev/null | head -n 1)
UNCLASS_COUNT=0
[ -n "$UNCLASS_BAM" ] && [ -f "$UNCLASS_BAM" ] && UNCLASS_COUNT=$(samtools view -c "$UNCLASS_BAM")

echo ""
echo "  Assigned (bc01–24) : $TOTAL_USED reads"
echo "  Unclassified       : $UNCLASS_COUNT reads"