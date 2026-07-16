#!/bin/bash
set -euo pipefail

SETUP_DIR="/scratch/users/k25053715/MSc_research/2_biopython_project/1_env_setup"
ENV_DIR="/scratch/users/k25053715/MSc_research/0_project_env/biopython_env"

mkdir -p "$SETUP_DIR"
mkdir -p "$ENV_DIR"

cd "$SETUP_DIR"

if command -v module >/dev/null 2>&1; then
    module load miniconda3 || module load anaconda3 || module load conda || true
fi

if ! command -v conda >/dev/null 2>&1; then
    echo "conda not found. Please load the correct conda module first."
    exit 1
fi

source "$(conda info --base)/etc/profile.d/conda.sh"

conda create -y --prefix "$ENV_DIR" python=3.11 biopython numpy
conda activate "$ENV_DIR"

python -c "import Bio, numpy; print('Biopython:', Bio.__version__)"
