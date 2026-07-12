import os
import openmm
import openmm.app as app
import openmm.unit as unit
from openmm import XmlSerializer

print("=== FGFR1 membrane minimization: COARSE stage ===")

INPUT_PDB = os.environ["INPUT_PDB"]
INPUT_XML = os.environ["SYSTEM_XML"]
OUTDIR = os.environ["OUTDIR"]

print(f"Input PDB: {INPUT_PDB}")
print(f"Input XML: {INPUT_XML}")
print(f"Output dir: {OUTDIR}")

# Step 1: load PDB + System XML
pdb = app.PDBFile(INPUT_PDB)
modeller = app.Modeller(pdb.topology, pdb.positions)
print(f"Loaded system: {modeller.topology.getNumAtoms()} atoms")

with open(INPUT_XML) as f:
    system = XmlSerializer.deserialize(f.read())
print(f"Loaded System from {INPUT_XML}")

# Step 2: add strong restraints on protein heavy atoms
restraint = openmm.CustomExternalForce(
    "0.5*k*((x-x0)^2 + (y-y0)^2 + (z-z0)^2)"
)
restraint.addGlobalParameter("k", 1000.0)
restraint.addPerParticleParameter("x0")
restraint.addPerParticleParameter("y0")
restraint.addPerParticleParameter("z0")

positions = modeller.positions
for atom in modeller.topology.atoms():
    if atom.residue.name not in ('POPC', 'POPE', 'DPPC', 'HOH', 'WAT', 'TIP3', 'TIP3P', 'Na+', 'Cl-'):
        if atom.element.symbol != app.element.hydrogen.symbol:
            idx = atom.index
            pos = positions[idx]
            restraint.addParticle(idx, pos)

system.addForce(restraint)

integrator = openmm.LangevinMiddleIntegrator(
    300 * unit.kelvin,
    1 / unit.picosecond,
    0.002 * unit.picoseconds
)

simulation = app.Simulation(modeller.topology, system, integrator)
simulation.context.setPositions(modeller.positions)

# Step 3: staged minimisation
state = simulation.context.getState(getEnergy=True)
energy = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
print(f"Initial energy: {energy:.1f} kJ/mol")

print("First pass: relax system (k=1000, maxIterations=2000, tol=1000 kJ/mol/nm)...")
simulation.minimizeEnergy(
    maxIterations=2000,
    tolerance=1000 * unit.kilojoule_per_mole / unit.nanometer
)

state = simulation.context.getState(getEnergy=True)
energy = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
print(f"After first pass: {energy:.1f} kJ/mol")

print("Second pass: tight minimization (k=100, maxIterations=1000, tol=10 kJ/mol/nm)...")
simulation.context.setParameter("k", 100.0)
simulation.minimizeEnergy(
    maxIterations=1000,
    tolerance=10 * unit.kilojoule_per_mole / unit.nanometer
)

state = simulation.context.getState(getEnergy=True, getPositions=True)
energy = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
positions = state.getPositions()
print(f"After second pass (coarse script final): {energy:.1f} kJ/mol")

# Step 4: save coarse-minimized structure
out_pdb = os.path.join(OUTDIR, '5_fgfr1_membrane_min_coarse_full.pdb')
app.PDBFile.writeFile(
    simulation.topology,
    positions,
    open(out_pdb, 'w')
)
print(f"Saved coarse-minimized full system to {out_pdb}")
