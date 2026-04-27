# Changelog

All notable changes to this project will be documented here.

## [Unreleased]

## [1.0.0] - 2026-04-27

### Added
- `00_setup_env.sh` — conda environment setup with minimap2, samtools, modkit
- `01_basecall_demux.sh` — Dorado SUP + 5mCG_5hmCG basecalling and SQK-RBK114-96 demultiplexing
- `02_align.sh` — Dorado aligner alignment to buffered reference FASTA; preserves MM/ML methylation tags
- `03_modkit_pileup.sh` — modkit pileup with `--cpg`, `--combine-strands`, `--modified-bases C:m`
- `04_build_beta_matrix.py` — Aggregates per-barcode BEDs into N_CpG × N_barcode beta matrix
- `run_all.sh` — Master orchestrator with `--from STEP` and `--min-cov` flags
- `README.md`, `environment.yml`, `.gitignore`

### Fixed
- Switched alignment from minimap2 to `dorado aligner` to prevent MM/ML tag loss
- Added `--no-trim` to dorado basecalling to preserve adapter sequences needed for demux
- Added `--combine-strands` to modkit to handle dual `_F`/`_R` contigs in buffered reference
