#!/bin/bash
#SBATCH --job-name=boltz2_install_gpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --gres=gpu:1
#SBATCH --partition=gpu
#SBATCH --time=00:45:00

set -euo pipefail

module purge
module load python/3.11.6-gcc-13.2.0 cuda/12.2.1-gcc-13.2.0

WORKDIR="/scratch/users/k25053715/MSc_research/0_project_env"
mkdir -p "$WORKDIR"
cd "$WORKDIR"

echo "=== Step 0: Clear cache ==="
rm -rf ~/.cache/torch ~/.boltz ~/.cache/huggingface

echo "=== Step 1: Fresh venv ==="
rm -rf boltz2_venv
python -m venv boltz2_venv
source boltz2_venv/bin/activate
pip install --upgrade pip setuptools wheel

echo "=== Step 2: Clone + install Boltz ==="
rm -rf boltz2
git clone https://github.com/jwohlwend/boltz.git boltz2
cd boltz2
pip install -e ".[cuda]"

echo "=== Step 3: GPU check (after install) ==="
nvidia-smi
python -c "import torch; print('CUDA:', torch.cuda.is_available()); print('Torch version:', torch.__version__); print('CUDA version:', torch.version.cuda)"

echo "=== Step 4: Test MSA + ligand (GPU) ==="
mkdir -p test_run_gpu
cat > test_run_gpu/pl.yaml << 'EOF'
version: 1
sequences:
  - protein:
      id: A
      sequence: "MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTQEDSYRKQVDSLLR"
  - ligand:
      id: B
      smiles: "CCO"
EOF

export TORCH_HOME="$WORKDIR/.cache/torch"
export HF_HOME="$WORKDIR/.cache/huggingface"
export XDG_CACHE_HOME="$WORKDIR/.cache"
srun boltz predict test_run_gpu/pl.yaml \
  --out_dir test_run_gpu/ \
  --use_msa_server \
  --accelerator gpu --devices 1

echo "GPU PIPELINE OK!"
ls -la test_run_gpu/predictions/ || echo "Check test_run_gpu/predictions/ for .cif files"

echo "=== USAGE (GPU job) ==="
echo "source $WORKDIR/boltz2_venv/bin/activate"
echo "boltz predict your_yaml.yaml --accelerator gpu --devices 1"
echo "DONE!"
