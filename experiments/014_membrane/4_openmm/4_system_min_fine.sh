#!/bin/bash
#SBATCH --job-name=3b_system_min_fine
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --partition=interruptible_cpu
#SBATCH --time=24:00:00
#SBATCH --requeue

module purge && module load anaconda3/2022.10-gcc-13.2.0
eval "$(conda shell.bash hook)"
conda activate /scratch/users/k25053715/MSc_research/0_project_env/openmm

INPUT_PDB="/scratch/users/k25053715/MSc_research/4_openmm_project/2_fgfr2_minimisation/014_membrane/014_5_fgfr2_membrane_min_coarse_full.pdb"

cp "$INPUT_PDB" ./014_5_fgfr2_membrane_min_coarse_full.pdb
echo "Copied: $INPUT_PDB → ./014_5_fgfr2_membrane_min_coarse_full.pdb"

cat > 4_system_min_fine.py << 'PYEOF'
import openmm
import openmm.app
import openmm.unit as unit
from openmm import XmlSerializer
import numpy as np

print("=== FGFR2 membrane minimization: FINE stage ===")

INPUT_PDB = '014_5_fgfr2_membrane_min_coarse_full.pdb'

# Step 1: load coarse-minimized PDB + same System XML
pdb = openmm.app.PDBFile(INPUT_PDB)
modeller = openmm.app.Modeller(pdb.topology, pdb.positions)
print(f"Loaded coarse-minimized system: {modeller.topology.getNumAtoms()} atoms")

with open('014_4_with_membrane_system.xml') as f:
    system = XmlSerializer.deserialize(f.read())
print("Loaded System from 014_4_with_membrane_system.xml")

# Step 2: weaker restraints on protein heavy atoms
restraint = openmm.CustomExternalForce(
    "0.5*k*((x-x0)^2 + (y-y0)^2 + (z-z0)^2)"
)
restraint.addGlobalParameter("k", 100.0)
restraint.addPerParticleParameter("x0")
restraint.addPerParticleParameter("y0")
restraint.addPerParticleParameter("z0")

positions = modeller.positions
for atom in modeller.topology.atoms():
    if atom.residue.name not in ('POPC', 'POPE', 'DPPC', 'HOH', 'WAT', 'TIP3', 'TIP3P', 'Na+', 'Cl-'):
        if atom.element.symbol != openmm.app.element.hydrogen:
            idx = atom.index
            pos = positions[idx]
            restraint.addParticle(idx, pos)

system.addForce(restraint)

integrator = openmm.LangevinMiddleIntegrator(
    300 * unit.kelvin,
    1 / unit.picosecond,
    0.002 * unit.picoseconds
)

simulation = openmm.app.Simulation(modeller.topology, system, integrator)
simulation.context.setPositions(modeller.positions)

# Step 3: fine minimisation
state = simulation.context.getState(getEnergy=True)
energy = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
print(f"Energy before fine: {energy:.1f} kJ/mol")

print("Fine min with weaker restraints (k=100, maxIterations=1000)...")
simulation.minimizeEnergy(
    tolerance=10 * unit.kilojoule_per_mole / unit.nanometer,
    maxIterations=1000
)

state = simulation.context.getState(getEnergy=True, getPositions=True)
energy = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
positions = state.getPositions()
print(f"Final energy after fine: {energy:.1f} kJ/mol")

# Step 4: save final minimized structures
openmm.app.PDBFile.writeFile(
    simulation.topology,
    positions,
    open('014_6_fgfr2_membrane_minimized_full.pdb', 'w')
)
print("Saved minimized full system to 014_6_fgfr2_membrane_minimized_full.pdb")

slim_min = openmm.app.Modeller(simulation.topology, positions)
waters_ions_min = [res for res in slim_min.topology.residues()
                   if res.name in ('HOH', 'WAT', 'TIP3', 'TIP3P', 'Na+', 'Cl-')]
slim_min.delete(waters_ions_min)
openmm.app.PDBFile.writeFile(
    slim_min.topology, slim_min.positions, open('014_6_fgfr2_membrane_minimized_no_water.pdb', 'w')
)
print("Saved minimized protein+membrane to 014_6_fgfr2_membrane_minimized_no_water.pdb")

PYEOF

echo "Starting FINE minimization from coarse-minimized membrane..."
python 4_system_min_fine.py

echo "Results (fine):"
ls -lh 014_6_fgfr2_membrane_minimized_*.pdb
