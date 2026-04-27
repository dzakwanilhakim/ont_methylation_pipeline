#!/bin/bash
# ============================================================
# 00_setup_env.sh
# Create and configure conda environment for ONT methylation pipeline
# ============================================================
set -euo pipefail

echo "[$(date)] Creating conda environment: ontmethyl_env"

conda create -n ontmethyl_env python=3.11 -y
conda activate ontmethyl_env

# minimap2
conda install -c bioconda -c conda-forge minimap2 -y
conda update minimap2 -y

# samtools
conda install -c bioconda -c conda-forge samtools -y
conda update samtools -y

# modkit
conda install -c conda-forge rust -y
cargo install --git https://github.com/nanoporetech/modkit.git
ln -s ~/.cargo/bin/modkit $(conda info --base)/envs/ontmethyl_env/bin/

# ncurses
conda install -c conda-forge ncurses --force-reinstall -y
conda update ncurses -y

# basic lib
conda install pandas numpy matplotlib seaborn -y


echo ""
echo "========================================"
echo " Environment setup complete."
echo " NOTE: Install Dorado manually (not on conda): https://software-docs.nanoporetech.com/dorado/latest/"
echo "   curl "https://cdn.oxfordnanoportal.com/software/analysis/dorado-1.4.0-linux-x64.tar.gz" -o dorado-1.4.0-linux-x64.tar.gz"
echo "   tar -xzf dorado-1.4.0-linux-x64.tar.gz"
echo "   dorado-1.4.0-linux-x64/bin/dorado --version"
echo "   dorado location :  "
echo "========================================"
