#!/bin/bash
#SBATCH --job-name=fgfr2_min_clean
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --partition=interruptible_cpu
#SBATCH --time=24:00:00

module purge && module load anaconda3/2022.10-gcc-13.2.0
eval "$(conda shell.bash hook)"
conda activate /scratch/users/k25053715/MSc_research/0_project_env/openmm

INPUT_PDB="/scratch/users/k25053715/MSc_research/2_biopython_project/2_straighten/012_finetune/012_FGFR2_clean_linkers.pdb"

cp "$INPUT_PDB" ./012_FGFR2_clean_linkers.pdb
echo "Copied: $INPUT_PDB → ./012_FGFR2_clean_linkers.pdb"

cat > fgfr2_minimize_clean.py << 'PYEOF'
import openmm
import openmm.app
import openmm.unit as unit
import pdbfixer

print("=== FGFR2 Clean Linkers Minimization ===")

############################################
# Step 1: Load and repair
############################################

fixer = pdbfixer.PDBFixer(filename='012_FGFR2_clean_linkers.pdb')
fixer.removeHeterogens(False)
fixer.findMissingResidues()
fixer.findMissingAtoms()
fixer.addMissingAtoms()

print(f"Repaired: {fixer.topology.getNumResidues()} residues")

############################################
# Save repaired
############################################

openmm.app.PDBFile.writeFile(fixer.topology, fixer.positions, open('012_1_clean_fixed.pdb', 'w'))

############################################
# Step 2: Forcefield + hydrogens
############################################

pdb = openmm.app.PDBFile('012_1_clean_fixed.pdb')
modeller = openmm.app.Modeller(pdb.topology, pdb.positions)

forcefield = openmm.app.ForceField('amber14-all.xml', 'implicit/obc2.xml')

disulfide_count = 0
for bond in list(modeller.topology.bonds()):
    a1, a2 = bond
    if a1.name == 'SG' and a2.name == 'SG':
        modeller.topology._bonds.remove(bond)
        disulfide_count += 1

print(f"Removed {disulfide_count} bad disulfides")

modeller.addHydrogens(forcefield, pH=7.4)
openmm.app.PDBFile.writeFile(modeller.topology, modeller.positions, open('012_2_clean_withH.pdb', 'w'))

############################################
# Step 3: System + simulation
############################################

system = forcefield.createSystem(
    modeller.topology,
    nonbondedMethod=openmm.app.CutoffNonPeriodic,
    nonbondedCutoff=1.0*unit.nanometer,
    constraints=openmm.app.HBonds,
    removeCMMotion=True
)

integrator = openmm.LangevinMiddleIntegrator(
    300*unit.kelvin,
    1/unit.picosecond,
    0.002*unit.picoseconds
)

simulation = openmm.app.Simulation(modeller.topology, system, integrator)
simulation.context.setPositions(modeller.positions)

state = simulation.context.getState(getEnergy=True)
energy = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
print(f"Initial energy: {energy:.1f} kJ/mol")

############################################
# Step 4: Progressive minimization
############################################

print("Coarse min...")
simulation.minimizeEnergy(
    tolerance=100*unit.kilojoule_per_mole/unit.nanometer,
    maxIterations=500
)

state = simulation.context.getState(getEnergy=True)
energy = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
print(f"After coarse: {energy:.1f} kJ/mol")

print("Fine min...")
simulation.minimizeEnergy(
    tolerance=1*unit.kilojoule_per_mole/unit.nanometer,
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
openmm.app.PDBFile.writeFile(
    simulation.topology,
    positions,
    open('012_3_fgfr2_clean_minimized.pdb', 'w')
)

final_energy = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
print(f"Saved minimized structure, final energy = {final_energy:.1f} kJ/mol")
PYEOF

echo "Starting minimization..."
python fgfr2_minimize_clean.py

echo "Results:"
ls -lh 012_*.pdb 012_FGFR2_clean_linkers.pdb
echo "Copy back: cp 012_3_fgfr2_clean_minimized.pdb /your/working/dir/"
