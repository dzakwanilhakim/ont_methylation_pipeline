# ONT Adaptive Sampling Methylation Pipeline

A complete end-to-end pipeline for **5mCG methylation analysis** using Oxford Nanopore Technology (ONT) adaptive sampling sequencing. Converts raw POD5 files to a beta methylation matrix across targeted CpG loci.

---

## Overview

```
POD5 (raw reads)
      │
      ▼
 01_basecall_demux.sh   — Dorado basecalling + demultiplexing (SQK-RBK114-96)
      │
      ▼
 02_align.sh            — Dorado aligner to buffered reference FASTA (preserves MM/ML tags)
      │
      ▼
 03_modkit_pileup.sh    — modkit pileup: per-barcode CpG methylation BEDs
      │
      ▼
 04_build_beta_matrix.py — Aggregate BEDs → beta matrix (N_CpGs × N_barcodes)
      │
      ▼
 beta_matrix_full.csv   — Final output
```

**Key design decisions:**
- Uses **Dorado aligner** (not minimap2) to preserve `MM`/`ML` methylation tags through alignment
- Uses a **buffered reference FASTA** (`5000 bp flanks` around each CpG target) for targeted adaptive sampling
- `--combine-strands` in modkit merges `_F` and `_R` contig hits into a single beta value per CpG

---

## Requirements

### Software

| Tool | Version | Install |
|------|---------|---------|
| [Dorado](https://github.com/nanoporetech/dorado) | ≥ 1.4.0 | Manual (see below) |
| [modkit](https://github.com/nanoporetech/modkit) | latest | via cargo |
| minimap2 | ≥ 2.24 | conda |
| samtools | ≥ 1.17 | conda |
| Python | ≥ 3.11 | conda |
| pandas, numpy | latest | conda/pip |

### Conda environment

```bash
bash src/00_setup_env.sh
conda activate ontmethyl_env
```

### Dorado (manual install)

Dorado is **not** on conda. Download and extract manually:

```bash
curl "https://cdn.oxfordnanoportal.com/software/analysis/dorado-1.4.0-linux-x64.tar.gz" \
    -o dorado-1.4.0-linux-x64.tar.gz
tar -xzf dorado-1.4.0-linux-x64.tar.gz
# Add to PATH or update DORADO= variable in src/01_basecall_demux.sh
```

---

## Project Structure

```
myproject/
├── data/
│   ├── pod5/                          ← raw POD5 input files
│   └── Buffered_5000_Fasta.fasta      ← buffered target reference (124 contigs for 62 CpGs)
├── analysis/                          ← all outputs (auto-created)
│   ├── 01_basecalled/                 ← calls.bam
│   ├── 02_demuxed/                    ← per-barcode BAMs
│   ├── 03_aligned/                    ← barcode01–24.sorted.bam
│   ├── 04_pileup/                     ← barcode01–24.bed
│   └── 05_beta_matrix/                ← final CSV outputs
└── src/
    ├── 00_setup_env.sh
    ├── 01_basecall_demux.sh
    ├── 02_align.sh
    ├── 03_modkit_pileup.sh
    ├── 04_build_beta_matrix.py
    └── run_all.sh
```

---

## Quick Start

### Run the full pipeline

```bash
bash src/run_all.sh
```

### Resume from a specific step

```bash
bash src/run_all.sh --from 2          # skip basecalling, start at alignment
bash src/run_all.sh --from 3          # skip to modkit pileup
bash src/run_all.sh --from 4 --min-cov 10   # rebuild beta matrix with stricter coverage
```

### Run individual steps

```bash
bash src/01_basecall_demux.sh
bash src/02_align.sh
bash src/03_modkit_pileup.sh
python src/04_build_beta_matrix.py --project /path/to/myproject --min-cov 5
```

---

## Configuration

Edit the variables at the top of each script before running:

### `01_basecall_demux.sh`
| Variable | Default | Description |
|----------|---------|-------------|
| `DORADO` | `dorado` | Path to dorado binary |
| `DORADO_MODEL` | `sup,5mCG_5hmCG` | Basecall + methylation model |
| `KIT` | `SQK-RBK114-96` | Barcoding kit |
| `FIRST_BC` / `LAST_BC` | `1` / `24` | Barcode range |

### `02_align.sh`
| Variable | Default | Description |
|----------|---------|-------------|
| `DEMUX_DIR` | (absolute path) | Path to demuxed BAMs |
| `THREADS` | `24` | CPU threads |

### `04_build_beta_matrix.py`
| Argument | Default | Description |
|----------|---------|-------------|
| `--min-cov` | `5` | Minimum valid reads to call a CpG |
| `--target-pos` | `5000` | 0-based CpG position in buffered FASTA |
| `--first-bc` / `--last-bc` | `1` / `24` | Barcode range |

---

## Outputs

All outputs are written to `analysis/05_beta_matrix/`:

| File | Description |
|------|-------------|
| `beta_matrix_full.csv` | N_CpG × 24 matrix; `NaN` = below coverage threshold |
| `beta_matrix_complete.csv` | CpGs with valid calls in **all** 24 barcodes |
| `beta_matrix_long.csv` | Tidy long format: `CpG_ID, barcode, beta, Nvalid` |
| `coverage_matrix.csv` | Read depth (Nvalid) per CpG per barcode |

---

## Reference FASTA Design

The buffered reference (`Buffered_5000_Fasta.fasta`) has **124 contigs** for **62 CpG targets**:

- Each CpG target generates **2 contigs**: `CG00001_F` and `CG00001_R`
- Each contig is `10,002 bp` total: `5000 bp flank + CpG + 5000 bp flank`
- The true CpG site is at **0-based position 5000** in each contig
- `--combine-strands` in modkit merges counts from `_F` and `_R` into one beta value

---

## Notes on MM/ML Tag Preservation

> **Important:** Use `dorado aligner` — not `minimap2` — for alignment.

`minimap2` converts BAM to SAM/FASTQ internally and **drops `MM`/`ML` methylation tags**. `dorado aligner` processes the BAM natively and preserves them. This is critical for `modkit pileup` to detect methylation calls.

---

## Citation / Context

This pipeline was developed for a Master's thesis project on HNSCC (Head and Neck Squamous Cell Carcinoma) DNA methylation classification using ONT adaptive sampling sequencing at Bandung Institute of Technology (ITB).

- Sequencing kit: SQK-RBK114-96 (24 barcoded patient samples)
- Target loci: 62 CpG sites selected from TCGA HNSCC methylation data
- Platform: Oxford Nanopore MinION (R10.4.1 flow cell)
- Downstream: Cross-platform ML classification (TCGA Illumina array → ONT validation)

---

## License

MIT
