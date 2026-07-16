#!/bin/bash
#SBATCH --job-name=1_protein_min
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --partition=interruptible_cpu
#SBATCH --time=24:00:00
#SBATCH --output=/scratch/users/k25053715/MSc_research/nextflow/2_pipeline/015_fgfr1/results/05_openmm/1_protein_min_%j.out
#SBATCH --error=/scratch/users/k25053715/MSc_research/nextflow/2_pipeline/015_fgfr1/results/05_openmm/1_protein_min_%j.err

set -euo pipefail

module purge
module load anaconda3/2022.10-gcc-13.2.0
eval "$(conda shell.bash hook)"
conda activate /scratch/users/k25053715/MSc_research/0_project_env/openmm

OUTDIR="/scratch/users/k25053715/MSc_research/nextflow/2_pipeline/015_fgfr1/results/05_openmm"
PYDIR="/scratch/users/k25053715/MSc_research/nextflow/2_pipeline/015_fgfr1/results/05_openmm/python"

mkdir -p "$OUTDIR" "$PYDIR"

start_ts=$(date +%s)
echo "[1_protein_min] Started at $(date)"

if [ "$#" -lt 1 ]; then
    echo "Usage: sbatch $0 /full/path/to/input.pdb"
    exit 1
fi

INPUT_PDB="$1"
export INPUT_PDB
export OUTDIR

if [ ! -f "$INPUT_PDB" ]; then
    echo "ERROR: Input PDB not found: $INPUT_PDB"
    exit 1
fi

echo "Using input PDB: $INPUT_PDB"
echo "Output directory: $OUTDIR"
echo "Python directory: $PYDIR"

cat > "$PYDIR/1_protein_min.py" << 'PYEOF'
import os
import openmm
import openmm.app as app
import openmm.unit as unit
import pdbfixer

INPUT_PDB = os.environ["INPUT_PDB"]
OUTDIR = os.environ["OUTDIR"]

print("=== FGFR1 Clean Linkers Minimization ===")
print(f"Input file: {INPUT_PDB}")
print(f"Output dir: {OUTDIR}")

############################################
# Step 1: Load and repair
############################################

fixer = pdbfixer.PDBFixer(filename=INPUT_PDB)
fixer.removeHeterogens(False)
fixer.findMissingResidues()
fixer.findMissingAtoms()
fixer.addMissingAtoms()

print(f"Repaired: {fixer.topology.getNumResidues()} residues")

############################################
# Save repaired
############################################

fixed_pdb = os.path.join(OUTDIR, "1_clean_fixed.pdb")
with open(fixed_pdb, "w") as f:
    app.PDBFile.writeFile(fixer.topology, fixer.positions, f)

############################################
# Step 2: Forcefield + hydrogens
############################################

pdb = app.PDBFile(fixed_pdb)
modeller = app.Modeller(pdb.topology, pdb.positions)

forcefield = app.ForceField("amber14-all.xml", "implicit/obc2.xml")

disulfide_count = 0
for bond in list(modeller.topology.bonds()):
    a1, a2 = bond
    if a1.name == "SG" and a2.name == "SG":
        modeller.topology._bonds.remove(bond)
        disulfide_count += 1

print(f"Removed {disulfide_count} bad disulfides")

modeller.addHydrogens(forcefield, pH=7.4)

withH_pdb = os.path.join(OUTDIR, "2_clean_withH.pdb")
with open(withH_pdb, "w") as f:
    app.PDBFile.writeFile(modeller.topology, modeller.positions, f)

############################################
# Step 3: System + simulation
############################################

system = forcefield.createSystem(
    modeller.topology,
    nonbondedMethod=app.CutoffNonPeriodic,
    nonbondedCutoff=1.0 * unit.nanometer,
    constraints=app.HBonds,
    removeCMMotion=True
)

integrator = openmm.LangevinMiddleIntegrator(
    300 * unit.kelvin,
    1 / unit.picosecond,
    0.002 * unit.picoseconds
)

simulation = app.Simulation(modeller.topology, system, integrator)
simulation.context.setPositions(modeller.positions)

state = simulation.context.getState(getEnergy=True)
energy = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
print(f"Initial energy: {energy:.1f} kJ/mol")

############################################
# Step 4: Progressive minimization
############################################

print("Coarse min...")
simulation.minimizeEnergy(
    tolerance=100 * unit.kilojoule_per_mole / unit.nanometer,
    maxIterations=500
)

state = simulation.context.getState(getEnergy=True)
energy = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
print(f"After coarse: {energy:.1f} kJ/mol")

print("Fine min...")
simulation.minimizeEnergy(
    tolerance=1 * unit.kilojoule_per_mole / unit.nanometer,
    maxIterations=2000
)

state = simulation.context.getState(getEnergy=True)
energy = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
print(f"Final energy: {energy:.1f} kJ/mol")

############################################
# Step 5: Save final structure
############################################

state = simulation.context.getState(getPositions=True, getEnergy=True)
positions = state.getPositions()

final_pdb = os.path.join(OUTDIR, "3_fgfr1_clean_minimized.pdb")
with open(final_pdb, "w") as f:
    app.PDBFile.writeFile(simulation.topology, positions, f)

final_energy = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
print(f"Saved minimized structure: {final_pdb}")
print(f"Final energy = {final_energy:.1f} kJ/mol")
PYEOF

echo "Starting minimization..."
python "$PYDIR/1_protein_min.py"

end_ts=$(date +%s)
elapsed=$(( end_ts - start_ts ))
printf "[1_protein_min] Finished in %d h %d min %d s\n" \
$((elapsed/3600)) $(((elapsed%3600)/60)) $((elapsed%60))

echo "Results:"
ls -lh "$OUTDIR"/*.pdb
echo "Python script saved at: $PYDIR/1_protein_min.py"
