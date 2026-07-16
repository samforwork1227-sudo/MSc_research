#!/bin/bash
#SBATCH --job-name=chimerax_pull
#SBATCH --output=install_chimerax_%j.out
#SBATCH --error=install_chimerax_%j.err
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G

set -euo pipefail

ENV_DIR=/scratch/users/k25053715/MSc_research/0_project_env/chimerax_env
cd "$ENV_DIR"

export SINGULARITY_CACHEDIR="$ENV_DIR/.singularity_cache"
export SINGULARITY_TMPDIR="$ENV_DIR/.singularity_tmp"
mkdir -p "$SINGULARITY_CACHEDIR" "$SINGULARITY_TMPDIR"

IMAGE="chimerax.sif"
URI="docker://maduprey/chimerax:1.5"

echo "Pulling $URI -> $IMAGE"
singularity pull "$IMAGE" "$URI"

echo "Pull success."
