#!/bin/bash
#SBATCH --job-name=clean_openmm
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=00:30:00

module purge && module load anaconda3/2022.10-gcc-13.2.0
eval "$(conda shell.bash hook)"

# Remove the old environment
conda env remove -n openmm -y

# Create a new environment with OpenMM and PDBFixer
conda create -n openmm python=3.10 openmm pdbfixer -c conda-forge -y

# Test
conda activate openmm
python -c "import openmm; print('OpenMM', openmm.__version__)"
python -c "import pdbfixer; print('PDBFixer OK')"
echo " Environment rebuild successful！"
